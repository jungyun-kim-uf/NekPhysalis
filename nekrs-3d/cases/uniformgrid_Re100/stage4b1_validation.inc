// Stage 4B-1: verification of pressure-lift operator and cubic CV balances.
// Stage5A: distributed over 4 or 6 MPI ranks; uniform affine HEX (P3..P10,h=8/7).
// IMPORTANT: the enclosing cube includes the fictitious constrained solid;
// its force estimate is an independent DIAGNOSTIC, not yet validated particle drag.
#include <limits>
#include <iomanip>

static int s4b1Every = 10;
static double s4b1HalfWidth = 16.0/7.0;
static dlong s4b1PressureCalls = 0;
static bool s4b1RecordPressure = false;
static double s4b1PreShiftL2 = 0, s4b1LiftAL2 = 0, s4b1RhsL2 = 0;
static std::vector<int> s4b1GLLDegree;
static std::vector<int> s4b1PressureDirichlet;
static std::vector<dlong> s4b1CVElements;
static std::vector<dfloat> s4b1X, s4b1Y, s4b1Z;
static std::vector<dfloat> s4b1Jw;
static bool s4b1HasPreviousMoment = false;
static double s4b1PreviousMoment = 0, s4b1PreviousTime = 0;

static bool s4b1Near(double a,double b) {return std::abs(a-b)<1e-9;}

// Mesh nodes are flattened as i + nq*(j+nq*k).
// This follows NekRS core/kernels/gradientVolumeHex3D.okl indexing.
static double s4b1Derivative(const std::vector<dfloat> &u,
                            dlong fieldOffset,dlong base,int nq,
                            int i,int j,int k,int field,int axis,
                            const dfloat *D)
{
  const int row=(axis==0?i:(axis==1?j:k));
  double r=0;
  for(int l=0;l<nq;++l){
    const int ii=(axis==0?l:i),jj=(axis==1?l:j),kk=(axis==2?l:k);
    r+=double(D[row*nq+l])*double(u[base+ii+nq*(jj+nq*kk)+field*fieldOffset]);
  }
  return (7.0/4.0)*r; // element spacing h=8/7
}

static void s4b1Setup()
{
  auto *mesh=nrs->meshV;
  platform->par->extract("casedata","verificationinterval",s4b1Every);
  platform->par->extract("casedata","controlvolumehalfwidth",s4b1HalfWidth);
  if(s4b1Every<1 || !s4b1Near(s4b1HalfWidth,16.0/7.0))
    throw std::runtime_error("Stage4B-1 requires validation interval>0 and CV halfwidth 2");
  if(mesh->Nq<4 || mesh->Nq>11 || mesh->Np!=mesh->Nq*mesh->Nq*mesh->Nq)
    throw std::runtime_error("Stage5C verifier requires uniform P3..P10");

  auto xyz=mesh->xyzHost();
  s4b1X=std::get<0>(xyz); s4b1Y=std::get<1>(xyz); s4b1Z=std::get<2>(xyz);
  s4b1Jw.assign(mesh->Nlocal,0);
  mesh->o_Jw.copyTo(s4b1Jw.data(),s4b1Jw.size());

  s4b1GLLDegree.assign(mesh->Nlocal,1);
  ogsGatherScatter(s4b1GLLDegree.data(),ogsInt,ogsAdd,mesh->ogs);

  const int nq=mesh->Nq, np=mesh->Np;
  for(dlong e=0;e<mesh->Nelements;++e){
    double mn[3]={1e100,1e100,1e100},mx[3]={-1e100,-1e100,-1e100};
    for(int n=0;n<np;++n){
      const dlong a=e*np+n;
      const double r[3]={double(s4b1X[a]),double(s4b1Y[a]),double(s4b1Z[a])};
      for(int d=0;d<3;++d){
        mn[d]=std::min(mn[d],r[d]);mx[d]=std::max(mx[d],r[d]);
      }
    }
    for(int d=0;d<3;++d)
      if(!s4b1Near(mx[d]-mn[d],8.0/7.0))
        throw std::runtime_error("Stage4B-1 CV requires all HEX elements h=8/7");
    bool inside=true;
    for(int d=0;d<3;++d)
      if(mn[d]<-s4b1HalfWidth-1e-10 ||
         mx[d]> s4b1HalfWidth+1e-10) inside=false;
    if(inside)s4b1CVElements.push_back(e);
  }
  long long localCv=static_cast<long long>(s4b1CVElements.size()),globalCv=0;
  MPI_Allreduce(&localCv,&globalCv,1,MPI_LONG_LONG,MPI_SUM,platform->comm.mpiComm());
  if(globalCv!=64)
    throw std::runtime_error("Stage5A globally expected exactly 64 HEX inside [-16/7,16/7]^3");

  // Check physical GLL differentiation on a linear coordinate field.
  double err=0;
  if(!s4b1CVElements.empty()){
    const dlong e=s4b1CVElements.front(),base=e*np;
    for(int k=0;k<nq;++k)for(int j=0;j<nq;++j)for(int i=0;i<nq;++i){
      for(int axis=0;axis<3;++axis){
        const auto &v=(axis==0?s4b1X:(axis==1?s4b1Y:s4b1Z));
        err=std::max(err,std::abs(s4b1Derivative(v,0,base,nq,i,j,k,0,axis,mesh->D)-1.0));
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&err,1,MPI_DOUBLE,MPI_MAX,platform->comm.mpiComm());
  if(err>1e-9)throw std::runtime_error("Stage5A distributed gradient check failed");

  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4b1_pressure_lift.csv",std::ios::trunc);
    f<<"solve_call,time,pressure_rhs_free_l2,ALp_free_l2,"
     <<"shifted_rhs_b_minus_ALp_free_l2,post_residual_b_minus_Ap_free_l2,"
     <<"post_residual_relative,masked_pressure_max_error,pressure_solver_final_residual\n";
    std::ofstream g("stage4b1_cv_balance.csv",std::ios::trunc);
    g<<"step,time,cube_halfwidth,cv_volume,mass_flux_signed,"
     <<"mass_flux_abs,normalized_mass_imbalance,divergence_l2,"
     <<"momentum_z,dmomentum_z_dt,advective_z_flux,pressure_z_flux,"
     <<"viscous_z_flux,Fz_CV_fictitious_domain_candidate,"
     <<"Cd_CV_fictitious_domain_candidate,Fz_Lamb_candidate,"
     <<"Cd_Lamb_candidate,Fz_CV_minus_Lamb\n";
    std::printf("[STAGE4B1] 64 HEX cubic control volume and linear-gradient test: PASS (%.3e)\n",err);
  }
}

static double s4b1AssembledFreeL2(std::vector<dfloat> a)
{
  auto *mesh=nrs->meshV;
  // Mask BEFORE assembly as in elliptic::_solve().
  for(dlong i=0;i<mesh->Nlocal;++i)if(s4b1PressureDirichlet[i])a[i]=0;
  ogsGatherScatter(a.data(),ogsDfloat,ogsAdd,mesh->ogs);
  double sum=0;
  for(dlong i=0;i<mesh->Nlocal;++i)
    sum+=double(a[i])*double(a[i])/s4b1GLLDegree[i];
  MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI_DOUBLE,MPI_SUM,platform->comm.mpiComm());
  return std::sqrt(sum);
}

// Called by extra Stage4B-1 source hook inside pressure solve with the
// ACTUAL unassembled pressure RHS and the pressure operator coefficient.
// A(Lp) and A(p) are evaluated through the exact NekRS elliptic Ax method.
static void s4b1PressureAudit(double time,int stage,
                              const occa::memory &o_rhs,
                              const occa::memory &o_lambda0,
                              bool afterSolve)
{
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal;
  if(!afterSolve){
    ++s4b1PressureCalls;
    s4b1RecordPressure=(s4b1PressureCalls%s4b1Every==0 || nrs->lastStep);
    if(!s4b1RecordPressure)return;
    s4b1PressureDirichlet.assign(N,0);
    const auto nMasked=nrs->fluid->ellipticSolverP->Nmasked();
    const auto o_mask=nrs->fluid->ellipticSolverP->o_maskIds();
    std::vector<dlong> h_mask(nMasked);
    if(nMasked>0)o_mask.copyTo(h_mask.data(),nMasked);
    for(auto m:h_mask) {
      if(m<0||m>=N)throw std::runtime_error("Stage4B-1 pressure mask ID out of range");
      s4b1PressureDirichlet[m]=1;
    }
    long long localMasked=static_cast<long long>(nMasked),globalMasked=0;
    MPI_Allreduce(&localMasked,&globalMasked,1,MPI_LONG_LONG,MPI_SUM,
                  platform->comm.mpiComm());
    if(globalMasked<1)
      throw std::runtime_error("Stage5A distributed pressure elliptic mask missing");
  }
  if(!s4b1RecordPressure)return;

  auto o_test=platform->deviceMemoryPool.reserve<dfloat>(N);
  auto o_A=platform->deviceMemoryPool.reserve<dfloat>(N);
  std::vector<dfloat> b(N),a(N),q(N);
  o_rhs.copyTo(b.data(),N);

  // Pure lift before solve; full solved pressure afterward.
  // Both use the same coefficients and local stiffness operator as NekRS.
  if(!afterSolve) o_test.copyFrom(s4b_ptarget.data(),N);
  const occa::memory &o_q = afterSolve ? nrs->fluid->o_P : o_test;
  nrs->fluid->ellipticSolverP->Ax(o_lambda0,o_NULL,o_q,o_A);
  o_A.copyTo(a.data(),N);
  std::vector<dfloat> r(N);
  for(dlong i=0;i<N;++i)r[i]=dfloat(double(b[i])-double(a[i]));
  const double rnorm=s4b1AssembledFreeL2(std::move(r));

  if(!afterSolve){
    s4b1PreShiftL2=rnorm;
    s4b1LiftAL2=s4b1AssembledFreeL2(a);
    s4b1RhsL2=s4b1AssembledFreeL2(b);
    if(!(std::isfinite(s4b1PreShiftL2)&&std::isfinite(s4b1LiftAL2)))
      throw std::runtime_error("Stage4B-1 nonfinite lift operator");
    return;
  }

  // The solve has produced p; check the original UNMASKED equation on all
  // free nodes using the assembled b-A(p), not the target vs itself.
  double maxMask=0;
  std::vector<dfloat> pp(N);
  nrs->fluid->o_P.copyTo(pp.data(),N);
  for(dlong i=0;i<N;++i)if(s4b_pmask[i])
    maxMask=std::max(maxMask,std::abs(double(pp[i])-double(s4b_ptarget[i])));
  MPI_Allreduce(MPI_IN_PLACE,&maxMask,1,MPI_DOUBLE,MPI_MAX,
                platform->comm.mpiComm());
  const double relative=rnorm/std::max(s4b1RhsL2,1e-30);
  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4b1_pressure_lift.csv",std::ios::app);
    f<<std::setprecision(17)<<s4b1PressureCalls<<","<<time<<","
     <<s4b1RhsL2<<","<<s4b1LiftAL2<<","<<s4b1PreShiftL2<<","
     <<rnorm<<","<<relative<<","<<maxMask<<","
     <<nrs->fluid->ellipticSolverP->finalResidual()<<"\n";
    std::printf("[STAGE4B1 PRESSURE] solve=%lld rfree/rhs=%.4e ALp=%.4e cageErr=%.3e\n",
                (long long)s4b1PressureCalls,relative,s4b1LiftAL2,maxMask);
  }
  if(!std::isfinite(rnorm)||!std::isfinite(relative)||!std::isfinite(maxMask))
    throw std::runtime_error("Stage4B-1 nonfinite pressure residual");
}

static void s4b1ControlVolume(int step,double time)
{
  if(step<=0 || (step%s4b1Every!=0 && !nrs->lastStep))return;
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal,off=nrs->fieldOffset;
  const int nq=mesh->Nq,np=mesh->Np;
  std::vector<dfloat> u(3*size_t(off)),p(N);
  nrs->fluid->o_U.copyTo(u.data(),u.size());
  nrs->fluid->o_P.copyTo(p.data(),p.size());

  const double mu=s3Viscosity,rho=double(RHO);
  double vol=0,moment=0,divSq=0,mass=0,massAbs=0;
  double adv=0,press=0,visc=0;
  const dfloat *D=mesh->D,*w=mesh->gllw;
  for(const auto e:s4b1CVElements){
    const dlong base=e*np;
    for(int k=0;k<nq;++k)for(int j=0;j<nq;++j)for(int i=0;i<nq;++i){
      const dlong a=base+i+nq*(j+nq*k);
      const double jw=s4b1Jw[a];
      vol+=jw;
      moment+=jw*rho*u[a+2*off];
      const double div=
        s4b1Derivative(u,off,base,nq,i,j,k,0,0,D)+
        s4b1Derivative(u,off,base,nq,i,j,k,1,1,D)+
        s4b1Derivative(u,off,base,nq,i,j,k,2,2,D);
      divSq+=jw*div*div;
    }
    // Only element faces ON the enclosing cubic CV boundary are integrated.
    // Each face belongs to exactly one interior CV element.
    for(int axis=0;axis<3;++axis)for(int side=0;side<2;++side){
      const int fixed=(side? nq-1:0);
      const int sign=(side? 1:-1);
      const dlong corner=base+(axis==0?fixed:0)+
                              nq*(axis==1?fixed:0)+nq*nq*(axis==2?fixed:0);
      const double coord=(axis==0?s4b1X[corner]:(axis==1?s4b1Y[corner]:s4b1Z[corner]));
      if(!s4b1Near(coord,sign*s4b1HalfWidth))continue;
      for(int r=0;r<nq;++r)for(int s=0;s<nq;++s){
        const int i=(axis==0?fixed:r);
        const int j=(axis==1?fixed:(axis==0?r:s));
        const int k=(axis==2?fixed:s);
        const dlong a=base+i+nq*(j+nq*k);
        const double wn=(16.0/49.0)*w[r]*w[s]; // (h/2)^2, h=8/7
        const double un=sign*u[a+axis*off],uz=u[a+2*off];
        mass+=wn*un;
        massAbs+=wn*std::abs(un);
        adv+=wn*rho*uz*un;
        press+=wn*p[a]*(axis==2?sign:0);
        const double duz_dn=sign*s4b1Derivative(u,off,base,nq,i,j,k,2,axis,D);
        const double dun_dz=sign*s4b1Derivative(u,off,base,nq,i,j,k,axis,2,D);
        visc+=wn*mu*(duz_dn+dun_dz);
      }
    }
  }
  // Stage5A: CV elements are distributed across ranks. Reduce ALL physical
  // integrals before the temporal finite difference and CSV output.
  double localTerms[8]={vol,moment,divSq,mass,massAbs,adv,press,visc};
  double globalTerms[8]={0,0,0,0,0,0,0,0};
  MPI_Allreduce(localTerms,globalTerms,8,MPI_DOUBLE,MPI_SUM,
                platform->comm.mpiComm());
  vol=globalTerms[0];moment=globalTerms[1];divSq=globalTerms[2];
  mass=globalTerms[3];massAbs=globalTerms[4];
  adv=globalTerms[5];press=globalTerms[6];visc=globalTerms[7];
  if(std::abs(vol-(64.0*512.0/343.0))>1e-8)
    throw std::runtime_error("Stage5A global CV volume integral differs from 64*(8/7)^3");
  double dmDt=std::numeric_limits<double>::quiet_NaN();
  double force=std::numeric_limits<double>::quiet_NaN();
  if(s4b1HasPreviousMoment && time>s4b1PreviousTime){
    dmDt=(moment-s4b1PreviousMoment)/(time-s4b1PreviousTime);
    // Fluid-on-particle force is minus particle-on-fluid momentum RHS.
    force=-(dmDt+adv+press-visc);
  }
  s4b1PreviousMoment=moment;
  s4b1PreviousTime=time;
  s4b1HasPreviousMoment=true;

  const auto c=s4BroadcastCoeffs(0);
  const auto fl=s4::lambForce(c,mu,rho);
  const double area=s3::PI*double(A)*double(A);
  const double denom=.5*rho*UINF*UINF*area;
  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4b1_cv_balance.csv",std::ios::app);
    f<<std::setprecision(17)<<step<<","<<time<<","<<s4b1HalfWidth<<","
     <<vol<<","<<mass<<","<<massAbs<<","<<std::abs(mass)/std::max(massAbs,1e-30)<<","
     <<std::sqrt(divSq/vol)<<","<<moment<<","<<dmDt<<","
     <<adv<<","<<press<<","<<visc<<","<<force<<","<<force/denom<<","
     <<fl[2]<<","<<fl[2]/denom<<","<<force-fl[2]<<"\n";
    std::printf("[STAGE4B1 CV] step=%d massFlux=%.3e divL2=%.3e CdCV*=%.6g CdLamb*=%.6g\n",
                step,mass,std::sqrt(divSq/vol),force/denom,fl[2]/denom);
  }
  if(!std::isfinite(vol)||!std::isfinite(mass)||!std::isfinite(adv)||
     !std::isfinite(press)||!std::isfinite(visc)||
     !std::isfinite(fl[2]))
    throw std::runtime_error("Stage4B-1 nonfinite CV diagnostic");
}
