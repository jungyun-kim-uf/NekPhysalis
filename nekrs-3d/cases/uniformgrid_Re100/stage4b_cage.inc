// Stage 4B: actual nodal Dirichlet cage targets for patched NekRS elliptic solvers.
// Stage5C: fixed sphere (0,0,0), a=1, P3-P10; run time validation required.
// Pressure is an analytical mean-free nodal lift; U is constrained on
// selected velocity cage + deep rigid interior.  No volume penalization is used.
#include <cmath>
#include <fstream>
#include <iomanip>

static deviceMemory<int> s4b_o_VMask, s4b_o_PMask;
static deviceMemory<dfloat> s4b_o_V, s4b_o_P;
static std::vector<int> s4b_vmask, s4b_pmask;
static std::vector<dfloat> s4b_vtarget, s4b_ptarget;
static int s4bUpdateCount=0;
static double s4bLastPressureMean=0;

static void s4bSetup()
{
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal, off=nrs->fieldOffset;
  s4b_vmask.assign(N,0);
  s4b_pmask.assign(N,0);
  s4b_vtarget.assign(3*size_t(off),dfloat(0));
  s4b_ptarget.assign(N,dfloat(0));
  const auto &g=stage2Cages.at(0);
  dlong locV=0,locP=0;
  for(dlong i=0;i<N;++i){
    s4b_vmask[i]=(g.vActive[i]||g.rigid[i])?1:0;
    s4b_pmask[i]=g.pActive[i]?1:0;
    locV+=s4b_vmask[i]; locP+=s4b_pmask[i];
  }
  dlong globV=0,globP=0;
  MPI_Allreduce(&locV,&globV,1,MPI_LONG_LONG,MPI_SUM,platform->comm.mpiComm());
  MPI_Allreduce(&locP,&globP,1,MPI_LONG_LONG,MPI_SUM,platform->comm.mpiComm());
  if(globV<=0 || globP<=0)
    throw std::runtime_error("STAGE5A globally empty velocity or pressure mask");
  s4b_o_VMask.resize(N);
  s4b_o_PMask.resize(N);
  s4b_o_V.resize(3*size_t(off));
  s4b_o_P.resize(N);
  s4b_o_VMask.copyFrom(s4b_vmask,s4b_vmask.size());
  s4b_o_PMask.copyFrom(s4b_pmask,s4b_pmask.size());
  s4b_o_V.copyFrom(s4b_vtarget,s4b_vtarget.size());
  s4b_o_P.copyFrom(s4b_ptarget,s4b_ptarget.size());

  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4b_constraints.csv",std::ios::trunc);
    f<<"step,time,velocity_cage_max_abs_error,velocity_cage_rms_error,"
       "pressure_cage_max_abs_error,pressure_cage_rms_error,"
       "pressure_target_weighted_mean\n";
    std::printf("[STAGE5A] global element-local V masks=%lld, P masks=%lld\n",
                (long long)globV,(long long)globP);
    std::printf("[STAGE4B] using patched elliptic nodal constraints; no penalization\n");
  }
}

// Invoked from PATCHED fluidSolver_t::solve before pressure and velocity solves,
// after its applyDirichlet() has assigned boundary values.
static void s4bApplyTargets(double,int)
{
  const dlong N=nrs->meshV->Nlocal, off=nrs->fieldOffset;
  np3dS4BApplyPressure(N,s4b_o_PMask,s4b_o_P,nrs->fluid->o_P);
  np3dS4BApplyVelocity(N,off,s4b_o_VMask,s4b_o_V,
                       nrs->fluid->o_U,nrs->fluid->o_Ue);
}

static void s4bCheckAfterSolve(int step,double time)
{
  if(step<=0) return;
  const dlong N=nrs->meshV->Nlocal,off=nrs->fieldOffset;
  std::vector<dfloat> u(3*size_t(off),0),p(N,0);
  nrs->fluid->o_U.copyTo(u.data(),u.size());
  nrs->fluid->o_P.copyTo(p.data(),p.size());
  double maxu=0,maxp=0,sumu=0,sump=0,wu=0,wp=0;
  for(dlong i=0;i<N;++i){
    if(s4b_vmask[i]){
      for(int k=0;k<3;++k){
        const double d=double(u[i+k*off])-double(s4b_vtarget[i+k*off]);
        maxu=std::max(maxu,std::abs(d));
        sumu+=s4Weights[i]*d*d;
      }
      wu+=s4Weights[i];
    }
    if(s4b_pmask[i]){
      const double d=double(p[i])-double(s4b_ptarget[i]);
      maxp=std::max(maxp,std::abs(d));
      sump+=s4Weights[i]*d*d;
      wp+=s4Weights[i];
    }
  }
  double localSums[4]={sumu,sump,wu,wp}, globalSums[4]={0,0,0,0};
  double localMaxes[2]={maxu,maxp},globalMaxes[2]={0,0};
  MPI_Allreduce(localSums,globalSums,4,MPI_DOUBLE,MPI_SUM,platform->comm.mpiComm());
  MPI_Allreduce(localMaxes,globalMaxes,2,MPI_DOUBLE,MPI_MAX,platform->comm.mpiComm());
  sumu=globalSums[0];sump=globalSums[1];
  wu=globalSums[2];wp=globalSums[3];
  maxu=globalMaxes[0];maxp=globalMaxes[1];
  if(wu<=0||wp<=0)throw std::runtime_error("Empty GLOBAL Stage5A checking weights");
  if(platform->comm.mpiRank()==0){
    std::ofstream f("stage4b_constraints.csv",std::ios::app);
    f<<std::setprecision(17)<<step<<","<<time<<","
     <<maxu<<","<<std::sqrt(sumu/(3*wu))<<","
     <<maxp<<","<<std::sqrt(sump/wp)<<","
     <<s4bLastPressureMean<<"\n";
    if(step%10==0||nrs->lastStep)
      std::printf("[STAGE4B] step=%d Ucage maxerr=%.4e Pcage maxerr=%.4e\n",
                  step,maxu,maxp);
  }
  if(!std::isfinite(maxu)||!std::isfinite(maxp))
    throw std::runtime_error("Stage4B nonfinite constrained field");
}

// After this step's solve, update targets from this step's GL coefficients.
// The next step uses them (single-pass time-lagged coupling).
static void s4bUpdateTargets(int step,double time)
{
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal,off=nrs->fieldOffset;
  const auto coeff=s4BroadcastCoeffs(0); // collective rank-0 coefficients
  const auto &part=particles.at(0);
  const auto state=s3Particle(part);
  const auto xyz=mesh->xyzHost();
  const auto &x=std::get<0>(xyz),&y=std::get<1>(xyz),&z=std::get<2>(xyz);
  double localPSum=0,localPW=0;
  std::fill(s4b_ptarget.begin(),s4b_ptarget.end(),dfloat(0));
  std::fill(s4b_vtarget.begin(),s4b_vtarget.end(),dfloat(0));
  for(dlong i=0;i<N;++i){
    const auto &g=stage2Cages.at(0);
    if(!g.vActive[i]&&!g.pActive[i]&&!g.rigid[i])continue;
    if(g.rigid[i]){
      for(int k=0;k<3;++k)s4b_vtarget[i+k*off]=dfloat(state.v[k]);
    }
    if(g.vActive[i]||g.pActive[i]){
      const auto val=s4::evaluate(
         {double(x[i]),double(y[i]),double(z[i])},state,coeff,
         double(RHO),s3Viscosity);
      if(g.vActive[i])
        for(int k=0;k<3;++k)s4b_vtarget[i+k*off]=dfloat(val.u[k]);
      if(g.pActive[i]){
        s4b_ptarget[i]=dfloat(val.p);
        localPSum+=s4Weights[i]*val.p;
        localPW+=s4Weights[i];
      }
    }
  }
  double globalPSum=0,globalPW=0;
  MPI_Allreduce(&localPSum,&globalPSum,1,MPI_DOUBLE,MPI_SUM,platform->comm.mpiComm());
  MPI_Allreduce(&localPW,&globalPW,1,MPI_DOUBLE,MPI_SUM,platform->comm.mpiComm());
  if(globalPW<=0)throw std::runtime_error("Stage4B missing Pin weights");
  const double mean=globalPSum/globalPW;
  s4bLastPressureMean=0.0;
  double remainder=0;
  for(dlong i=0;i<N;++i) if(s4b_pmask[i]){
    s4b_ptarget[i]=dfloat(double(s4b_ptarget[i])-mean);
    remainder+=s4Weights[i]*s4b_ptarget[i];
  }
  MPI_Allreduce(MPI_IN_PLACE,&remainder,1,MPI_DOUBLE,MPI_SUM,platform->comm.mpiComm());
  s4bLastPressureMean=remainder/globalPW;
  if(std::abs(s4bLastPressureMean)>1e-8)
    throw std::runtime_error("Stage4B pressure target is not weighted mean-free");
  s4b_o_V.copyFrom(s4b_vtarget,s4b_vtarget.size());
  s4b_o_P.copyFrom(s4b_ptarget,s4b_ptarget.size());
  s4bUpdateCount++;
  if(step%10==0 && platform->comm.mpiRank()==0)
    printf("[STAGE4B] updated lagged targets at step %d; pressure mean=%.3e\n",
           step,s4bLastPressureMean);
}
