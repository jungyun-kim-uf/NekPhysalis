// Stage4B-1: secondary diagnostics of newly updated Lamb targets.
// The elliptic cage constraints are applied by Stage4B in the patched solver.
// Root GL coefficients are broadcast so every MPI rank computes its own local targets.
#include "stage4_reconstruct.hpp"
#include <iomanip>

static std::vector<dfloat> s4Weights;
static s3::V s4dummy;
static deviceMemory<dfloat> s4_o_PTarget; // prepared on GPU, diagnostics only; not coupled
static deviceMemory<dfloat> s4_o_UTarget; // prepared on GPU, diagnostics only; not coupled
static std::vector<dfloat> s4_ptarget;
static std::vector<dfloat> s4_utarget;
static int s4Every=10;

static std::vector<s3::Coeff> s4BroadcastCoeffs(size_t pid)
{
  const int rank=platform->comm.mpiRank();
  const auto nm=s3::modes().size();
  std::vector<double> packed(6*nm,0.0);
  if(rank==0){
    for(size_t j=0;j<nm;++j){
      const auto &c=s3Committed.at(pid).at(j);
      packed[6*j+0]=c.p.real();
      packed[6*j+1]=c.p.imag();
      packed[6*j+2]=c.phi.real();
      packed[6*j+3]=c.phi.imag();
      packed[6*j+4]=c.v.real();
      packed[6*j+5]=c.v.imag();
    }
  }
  MPI_Bcast(packed.data(),int(packed.size()),MPI_DOUBLE,0,platform->comm.mpiComm());
  std::vector<s3::Coeff> out(nm);
  for(size_t j=0;j<nm;++j){
    out[j].p={packed[6*j+0],packed[6*j+1]};
    out[j].phi={packed[6*j+2],packed[6*j+3]};
    out[j].v={packed[6*j+4],packed[6*j+5]};
  }
  return out;
}

static void s4Setup()
{
  if(samplingEvery!=1)
    throw std::runtime_error("Stage 4A requires glSamplingInterval=1 (single-pass coefficient updates)");
  auto *mesh=nrs->meshV;
  s4Weights.resize(mesh->Nlocal);
  mesh->o_Jw.copyTo(s4Weights.data(),s4Weights.size());
  const size_t velSize=3*size_t(nrs->fieldOffset);
  s4_ptarget.assign(mesh->Nlocal,0);
  s4_utarget.assign(velSize,0);
  s4_o_PTarget.resize(mesh->Nlocal);
  s4_o_UTarget.resize(velSize);
  const double ns=s4::noSlipManufacturedTest();
  const double recon=s4::reconstructionConsistencyTest();
  if(!(std::isfinite(ns) && ns<1.e-11 && std::isfinite(recon) && recon<1.e-10))
    throw std::runtime_error("Stage4 Lamb no-slip / reconstruction consistency test FAILED");
  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4_analytic_tests.csv",std::ios::trunc);
    f<<"test,error,threshold,pass\n";
    f<<std::setprecision(17)
     <<"manufactured_no_slip,"<<ns<<",1e-11,1\n"
     <<"pressure_tangential_surface_match,"<<recon<<",1e-10,1\n";
    std::ofstream s("stage4_cage_targets_summary.csv",std::ios::trunc);
    s<<"step,time,particle_id,"
       "max_U_cage_mismatch,rms_U_cage_mismatch,"
       "weighted_mean_raw_p_correction,weighted_mean_free_lift,"
       "max_abs_pressure_lift,Fz_Lamb_candidate_NOT_VALIDATED\n";
    std::printf("[STAGE4A] no-slip error %.4e; pressure/tangential consistency %.4e\n",ns,recon);
    std::printf("[STAGE4B1] Updated-target diagnostic; patched core constraints active.\n");
  }
}

static void s4Diagnostics(int step,double time)
{
  if(step%s4Every!=0 && !nrs->lastStep)return;
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal, off=nrs->fieldOffset;
  const int rank=platform->comm.mpiRank();
  const MPI_Comm comm=platform->comm.mpiComm();
  const auto xyz=mesh->xyzHost();
  const auto &x=std::get<0>(xyz), &y=std::get<1>(xyz), &z=std::get<2>(xyz);
  std::vector<dfloat> fluidP(N,0),fluidU(3*size_t(off),0);
  nrs->fluid->o_P.copyTo(fluidP.data(),fluidP.size());
  nrs->fluid->o_U.copyTo(fluidU.data(),fluidU.size());

  std::fill(s4_utarget.begin(),s4_utarget.end(),dfloat(0));
  std::fill(s4_ptarget.begin(),s4_ptarget.end(),dfloat(0));

  for(size_t pid=0;pid<particles.size();++pid){
    const auto c=s4BroadcastCoeffs(pid);
    const auto state=s3Particle(particles[pid]);
    const auto &cage=stage2Cages.at(pid);
    double locUMax=0,locUSq=0,locUCnt=0;
    double locPSum=0,locPW=0;
    std::vector<double> rawCorrection(N,0.0);
    for(dlong i=0;i<N;++i){
      if(!cage.vActive[i] && !cage.pActive[i])continue;
      const auto value=s4::evaluate({double(x[i]),double(y[i]),double(z[i])},
                                     state,c,double(RHO),s3Viscosity);
      if(cage.vActive[i]){
        for(int k=0;k<3;++k){
          const double err=double(fluidU[i+k*off])-value.u[k];
          locUMax=std::max(locUMax,std::abs(err));
          locUSq+=s4Weights[i]*err*err;
          s4_utarget[i+k*off]=dfloat(value.u[k]);
        }
        locUCnt+=s4Weights[i];
      }
      if(cage.pActive[i]){
        rawCorrection[i]=value.p-fluidP[i];
        locPSum+=s4Weights[i]*rawCorrection[i];
        locPW+=s4Weights[i];
      }
    }
    double globPSum=0,globPW=0,globUMax=0,globUSq=0,globUW=0;
    MPI_Allreduce(&locPSum,&globPSum,1,MPI_DOUBLE,MPI_SUM,comm);
    MPI_Allreduce(&locPW,&globPW,1,MPI_DOUBLE,MPI_SUM,comm);
    MPI_Allreduce(&locUMax,&globUMax,1,MPI_DOUBLE,MPI_MAX,comm);
    MPI_Allreduce(&locUSq,&globUSq,1,MPI_DOUBLE,MPI_SUM,comm);
    MPI_Allreduce(&locUCnt,&globUW,1,MPI_DOUBLE,MPI_SUM,comm);
    if(globPW<=0 || globUW<=0)
      throw std::runtime_error("Stage4A cage has zero mass weights");
    const double mean=globPSum/globPW;
    double locLiftSum=0,locLiftMax=0;
    for(dlong i=0;i<N;++i) if(cage.pActive[i]){
      const double lift=rawCorrection[i]-mean;
      s4_ptarget[i]=dfloat(lift);
      locLiftSum+=s4Weights[i]*lift;
      locLiftMax=std::max(locLiftMax,std::abs(lift));
    }
    double globalLiftSum=0,globalLiftMax=0;
    MPI_Allreduce(&locLiftSum,&globalLiftSum,1,MPI_DOUBLE,MPI_SUM,comm);
    MPI_Allreduce(&locLiftMax,&globalLiftMax,1,MPI_DOUBLE,MPI_MAX,comm);
    const double actualLiftMean=globalLiftSum/globPW;
    if(std::abs(actualLiftMean)>1.e-8*std::max(1.0,globalLiftMax))
      throw std::runtime_error("Stage4A mean-free pressure lift verification FAILED");

    const auto F=s4::lambForce(c,s3Viscosity,double(RHO));
    if(rank==0){
      std::ofstream f("stage4_cage_targets_summary.csv",std::ios::app);
      f<<std::setprecision(17)
       <<step<<","<<time<<","<<particles[pid].id<<","
       <<globUMax<<","
       <<std::sqrt(globUSq/(3*globUW))<<","
       <<mean<<","<<actualLiftMean<<","<<globalLiftMax<<","<<F[2]<<"\n";
      std::printf("[STAGE4A] step=%d particle=%d Ucage mismatch max=%.4e; "
                  "lift mean=%.4e; |lift|max=%.4e; Fz(candidate)=%.4e (candidate, NOT VALIDATED DRAG)\n",
                  step,particles[pid].id,globUMax,actualLiftMean,globalLiftMax,F[2]);
    }
    if(nrs->lastStep){
      char fname[128];
      std::snprintf(fname,sizeof(fname),"stage4_targets_rank%04d_particle%04d.csv",
                    rank,particles[pid].id);
      std::ofstream f(fname,std::ios::trunc);
      f<<std::setprecision(17);
      f<<"rank,local_node_index,x,y,z,Uout,Pin,"
         "ux_target,uy_target,uz_target,pressure_meanfree_lift\n";
      for(dlong i=0;i<N;++i)if(cage.vActive[i]||cage.pActive[i]){
        f<<rank<<","<<i<<","<<x[i]<<","<<y[i]<<","<<z[i]<<","
         <<cage.vActive[i]<<","<<cage.pActive[i]<<","
         <<s4_utarget[i]<<","<<s4_utarget[i+off]<<","<<s4_utarget[i+2*off]
         <<","<<s4_ptarget[i]<<"\n";
      }
    }
  }
  // Secondary diagnostic targets; Stage4B supplies actual solver constraints.
  s4_o_PTarget.copyFrom(s4_ptarget,s4_ptarget.size());
  s4_o_UTarget.copyFrom(s4_utarget,s4_utarget.size());
}
