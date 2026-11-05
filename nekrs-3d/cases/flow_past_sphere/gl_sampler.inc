// Stage 3: actual 3D (n,m) GL48x32 scalar products and p,phi,v Lamb
// coefficient extraction, using distributed NekRS findpts interpolation.
// Stage4B-1: targets are applied by patched elliptic masks; no volume penalty.
// This include is HOST ONLY; needs Particle, particles, RHO defined in np3d.udf.
#include <memory>
#include <iomanip>
#include "stage3_lamb_math.hpp"

static std::unique_ptr<pointInterpolation_t> s3Interpolator;
static std::vector<dfloat> s3X,s3Y,s3Z;
static std::vector<s3::Point> s3Points; // root only, concatenated per particle
static deviceMemory<dfloat> s3_o_U, s3_o_P;
static int samplingEvery=10;
static dfloat qOverA=1.14, coefficientBeta=0.01;
static std::vector<std::vector<s3::Coeff>> s3Committed; // per particle (root)
static double s3Viscosity=0;

static s3::ParticleState s3Particle(const Particle &p)
{
  return {p.id, {double(p.x),double(p.y),double(p.z)},
                {double(p.vx),double(p.vy),double(p.vz)},double(p.radius)};
}

static void glSetup()
{
  if (samplingEvery<=0)
    throw std::runtime_error("Stage3 requires glSamplingInterval>0");
  if (!(qOverA>1.0 && qOverA<=2.0 &&
        s3Viscosity>0 && coefficientBeta>0 && coefficientBeta<=1))
    throw std::runtime_error("Stage3 invalid q/a, viscosity or coefficient beta");
  if(platform->comm.mpiRank()==0) {
    s3Committed.resize(particles.size(), std::vector<s3::Coeff>(s3::modes().size()));
    for(const auto &part:particles) {
      auto qpoints=s3::spherePoints(s3Particle(part),qOverA);
      for(const auto &point:qpoints) {
        s3Points.push_back(point);
        s3X.push_back(point.x[0]);
        s3Y.push_back(point.x[1]);
        s3Z.push_back(point.x[2]);
      }
    }
    s3_o_U.resize(3*s3Points.size());
    s3_o_P.resize(s3Points.size());
    const double err=s3::manufacturedTest();
    std::printf("[STAGE3] GL48x32 synthetic (p,phi,v) roundtrip max error=%.5e\n",err);
    if(!(std::isfinite(err) && err<1e-10))
      throw std::runtime_error("Stage3 manufactured coefficient inversion FAILED");
    {
      std::ofstream f("stage3_manufactured_test.csv",std::ios::trunc);
      f << "q_over_a,N,gl_mu,gl_phi,max_coefficient_error,pass\n"
        << qOverA << "," << s3::LMAX << "," << s3::NMU << "," << s3::NPHI << "," << std::setprecision(17) << err << ",1\n";
    }
    {
      std::ofstream f("stage3_lamb_coefficients.csv",std::ios::trunc);
      f << "step,time,particle_id,n,m,"
           "p_raw_re,p_raw_im,phi_raw_re,phi_raw_im,v_raw_re,v_raw_im,"
           "p_relaxed_re,p_relaxed_im,phi_relaxed_re,phi_relaxed_im,"
           "v_relaxed_re,v_relaxed_im\n";
    }
    {
      std::ofstream f("stage3_sampling_summary.csv",std::ios::trunc);
      f << "step,time,particle_id,num_samples,mean_pressure,"
           "max_raw_coefficient_abs,max_relaxed_change_abs,"
           "point_search_failure_count\n";
    }
    std::printf("[STAGE3] root query points=%zu (GL48x32 per particle), q/a=%.6f, "
                "beta=%.4f, N=%d\n",s3Points.size(),double(qOverA),double(coefficientBeta),s3::LMAX);
  } else {
    // Stage5A MPI: keep a VALID point on non-root ranks so that OCCA
    // setPoints()/eval() never receives a zero-byte query/output buffer.
    // These rank-local dummy samples are discarded; only root's 1536
    // sphere samples enter the spherical-harmonic projection.
    s3X.push_back(dfloat(0));s3Y.push_back(dfloat(0));s3Z.push_back(dfloat(0));
    s3_o_U.resize(3);s3_o_P.resize(1);
  }

  s3Interpolator.reset(new pointInterpolation_t(nrs->meshV,platform->comm.mpiComm()));
  s3Interpolator->setPoints(s3X,s3Y,s3Z);
  s3Interpolator->find();

  int nfail=0;
  for(const auto code:s3Interpolator->data().code)
    if(code==pointInterpolation_t::CODE_NOT_FOUND)++nfail;
  int globalFail=0;
  MPI_Allreduce(&nfail,&globalFail,1,MPI_INT,MPI_SUM,platform->comm.mpiComm());
  if(globalFail!=0)
    throw std::runtime_error("Stage3 off-grid GL48x32 point search failed");
  if(platform->comm.mpiRank()==0)
    std::printf("[STAGE3] All %zu sampling points located successfully.\n",s3Points.size());
}

static void glSample(int step,double time)
{
  if(samplingEvery<=0 || step%samplingEvery!=0)return;
  // Interpolation evaluation is collective; ALL ranks must enter here.
  // Root owns physical GL points; other ranks carry one in-domain dummy.
  const dlong count=s3X.size();
  s3Interpolator->eval(3,nrs->fluid->fieldOffset,nrs->fluid->o_U,
                       count,s3_o_U);
  s3Interpolator->eval(1,nrs->fluid->fieldOffset,nrs->fluid->o_P,
                       count,s3_o_P);
  if(platform->comm.mpiRank()!=0)return;

  // User's verified NekRS v26.0.1 typed-deviceMemory API: vector, NOT .data().
  std::vector<dfloat> hp(count), hu(3*count);
  s3_o_P.copyTo(hp,count);
  s3_o_U.copyTo(hu,3*count);

  const auto ml=s3::modes();
  std::ofstream f("stage3_lamb_coefficients.csv",std::ios::app);
  std::ofstream summary("stage3_sampling_summary.csv",std::ios::app);
  f << std::setprecision(16);
  summary << std::setprecision(16);

  const size_t nPerParticle=s3::NMU*s3::NPHI;
  for(size_t pid=0;pid<particles.size();++pid) {
    const size_t begin=pid*nPerParticle;
    std::vector<s3::Point> qpoints(s3Points.begin()+begin,
                                   s3Points.begin()+begin+nPerParticle);
    std::vector<s3::Fields> samples(nPerParticle);
    for(size_t k=0;k<nPerParticle;++k){
      const size_t index=begin+k;
      samples[k].p=hp[index];
      for(int d=0;d<3;++d)samples[k].u[d]=hu[index+d*size_t(count)];
    }
    double pMean=0;
    const auto products=s3::project(qpoints,samples,s3Particle(particles[pid]),&pMean);
    const auto raw=s3::invert(products,s3Particle(particles[pid]),
                              double(qOverA),double(RHO),s3Viscosity);
    double maxRaw=0,maxDelta=0;
    for(size_t j=0;j<raw.size();++j) {
      auto &c=s3Committed[pid][j];
      const auto prev=c;
      c.p=(1-double(coefficientBeta))*c.p+double(coefficientBeta)*raw[j].p;
      c.phi=(1-double(coefficientBeta))*c.phi+double(coefficientBeta)*raw[j].phi;
      c.v=(1-double(coefficientBeta))*c.v+double(coefficientBeta)*raw[j].v;
      maxRaw=std::max({maxRaw,std::abs(raw[j].p),
                       std::abs(raw[j].phi),std::abs(raw[j].v)});
      maxDelta=std::max({maxDelta,std::abs(c.p-prev.p),
                         std::abs(c.phi-prev.phi),std::abs(c.v-prev.v)});
      f << step << ',' << time << ',' << particles[pid].id << ','
        << ml[j].n << ',' << ml[j].m << ','
        << raw[j].p.real() << ',' << raw[j].p.imag() << ','
        << raw[j].phi.real() << ',' << raw[j].phi.imag() << ','
        << raw[j].v.real() << ',' << raw[j].v.imag() << ','
        << c.p.real() << ',' << c.p.imag() << ','
        << c.phi.real() << ',' << c.phi.imag() << ','
        << c.v.real() << ',' << c.v.imag() << '\n';
    }
    summary << step << ',' << time << ',' << particles[pid].id << ','
            << nPerParticle << ',' << pMean << ',' << maxRaw
            << ',' << maxDelta << ",0\n";
    std::printf("[STAGE3] step=%d particle=%d pMean=%g max|raw coeff|=%g "
                "max|relax increment|=%g (cage-constrained Stage4B field; Cd not validated)\n",
                step,particles[pid].id,pMean,maxRaw,maxDelta);
  }
}
