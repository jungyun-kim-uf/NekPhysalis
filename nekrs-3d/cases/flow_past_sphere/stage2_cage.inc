// P6 h=8/7 Physalis cage topology sweep: S1/S2/S3 and C2/C3.
// Graph edges are tensor-neighbor GLL pairs, globally synchronized by OGS.
// This case-specific source is intended for the existing patched NekRS solver.
// Uin cage: analytic target at selected inner layers; deeper interior rigid.
// Uout cage: analytic target at selected outer layers; all interior rigid.
// Pressure cage: analytic target at selected inner or outer layers.
#include "ogs.hpp"
#include <array>
#include <limits>
#include <tuple>
#include <set>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <cctype>

struct Stage2Cage {
  int particleID = -1;
  std::vector<int> uoutS1, pinS1, interior;
  std::vector<int> vActive, pActive, rigid;
};
static std::vector<Stage2Cage> stage2Cages;

static int s2Parse(const std::string &text, bool isVel,
                   bool &isInside, std::array<int,3> &use)
{
  std::string t=text;
  std::transform(t.begin(),t.end(),t.begin(),
      [](unsigned char c){return char(std::tolower(c));});
  const std::string prefixIn=isVel?"uin-":"pin-";
  const std::string prefixOut=isVel?"uout-":"pout-";
  std::string suffix;
  if(t.rfind(prefixIn,0)==0){isInside=true;suffix=t.substr(prefixIn.size());}
  else if(t.rfind(prefixOut,0)==0){isInside=false;suffix=t.substr(prefixOut.size());}
  else throw std::runtime_error("Invalid cage side: "+text);
  use={0,0,0};
  if(suffix=="s1") use[0]=1;
  else if(suffix=="s2")use[1]=1;
  else if(suffix=="s3")use[2]=1;
  else if(suffix=="c2")use={1,1,0};
  else if(suffix=="c3")use={1,1,1};
  else throw std::runtime_error("Invalid cage layer "+text+" (S1/S2/S3/C2/C3)");
  return use[2] ? 3 : (use[1] ? 2 : 1);
}

static void cageStage2Setup()
{
  auto *mesh=nrs->meshV;
  const dlong N=mesh->Nlocal;
  const int nq=mesh->Nq,np=mesh->Np;
  const int rank=platform->comm.mpiRank();
  const MPI_Comm comm=platform->comm.mpiComm();
  if(mesh->dim!=3 || nq!=7 || np!=343)
    throw std::runtime_error("Topology sweep requires 3D P6 HEX");

  std::string vel="Uout-S1",prs="Pin-S1";
  int harmonicOrder=3,glTheta=32,glPhi=32;
  platform->par->extract("casedata","cagevelocity",vel);
  platform->par->extract("casedata","cagepressure",prs);
  platform->par->extract("casedata","harmonicorder",harmonicOrder);
  platform->par->extract("casedata","glthetapoints",glTheta);
  platform->par->extract("casedata","glphipoints",glPhi);
  if(harmonicOrder!=s3::LMAX||glTheta!=s3::NMU||glPhi!=s3::NPHI)
    throw std::runtime_error("N-sweep parameters disagree with compiled harmonic basis and GL grid");
  bool vInside=false,pInside=true;
  std::array<int,3> vUse{},pUse{};
  const int vDepth=s2Parse(vel,true,vInside,vUse);
  (void)s2Parse(prs,false,pInside,pUse);

  if(particles.size()!=1)throw std::runtime_error("Exactly one particle required");
  const auto &part=particles.at(0);
  const auto xyz=mesh->xyzHost();
  const auto &x=std::get<0>(xyz);
  const auto &y=std::get<1>(xyz);
  const auto &z=std::get<2>(xyz);
  Stage2Cage cage;
  cage.particleID=part.id;
  cage.interior.assign(N,0);
  cage.vActive.assign(N,0);
  cage.pActive.assign(N,0);
  cage.rigid.assign(N,0);
  std::array<std::vector<int>,3> outward,inward;
  for(int k=0;k<3;++k){outward[k].assign(N,0);inward[k].assign(N,0);}
  for(dlong i=0;i<N;++i){
    const double dx=x[i]-part.x,dy=y[i]-part.y,dz=z[i]-part.z;
    const double r=std::sqrt(dx*dx+dy*dy+dz*dz)/part.radius;
    cage.interior[i]=r<1.0-1e-12 ? 1:0;
  }
  // Shared-node classification must agree across all elements and MPI ranks.
  ogsGatherScatter(cage.interior.data(),ogsInt,ogsMax,mesh->ogs);

  // For each element edge, visit both directions.
  auto visitEdges=[&](auto fn){
    for(dlong e=0;e<mesh->Nelements;++e){
      const dlong base=e*np;
      for(int k=0;k<nq;++k)for(int j=0;j<nq;++j)for(int i=0;i<nq;++i){
        const dlong a=base+i+nq*(j+nq*k);
        if(i+1<nq)fn(a,a+1);
        if(j+1<nq)fn(a,a+nq);
        if(k+1<nq)fn(a,a+nq*nq);
      }
    }
  };
  visitEdges([&](dlong a,dlong b){
    if(cage.interior[a]!=cage.interior[b]){
      const dlong ia=cage.interior[a]?a:b;
      const dlong oa=cage.interior[a]?b:a;
      inward[0][ia]=1;outward[0][oa]=1;
    }
  });
  ogsGatherScatter(outward[0].data(),ogsInt,ogsMax,mesh->ogs);
  ogsGatherScatter(inward[0].data(),ogsInt,ogsMax,mesh->ogs);
  for(int k=1;k<3;++k){
    visitEdges([&](dlong a,dlong b){
      if(cage.interior[a]!=cage.interior[b])return;
      auto &curr=cage.interior[a]?inward:outward;
      if(curr[k-1][a])curr[k][b]=1;
      if(curr[k-1][b])curr[k][a]=1;
    });
    for(dlong i=0;i<N;++i){
      for(int earlier=0;earlier<k;++earlier){
        if(outward[earlier][i])outward[k][i]=0;
        if(inward[earlier][i])inward[k][i]=0;
      }
    }
    ogsGatherScatter(outward[k].data(),ogsInt,ogsMax,mesh->ogs);
    ogsGatherScatter(inward[k].data(),ogsInt,ogsMax,mesh->ogs);
  }
  cage.uoutS1=outward[0];
  cage.pinS1=inward[0];

  for(dlong i=0;i<N;++i) {
    const auto &vls=vInside?inward:outward;
    const auto &pls=pInside?inward:outward;
    for(int k=0;k<3;++k){
      if(vUse[k]&&vls[k][i])cage.vActive[i]=1;
      if(pUse[k]&&pls[k][i])cage.pActive[i]=1;
    }
    if(!vInside) {
      cage.rigid[i]=cage.interior[i];
    } else if(cage.interior[i]) {
      // Rigid only *behind* selected inner cage (deeper than its maximum layer).
      bool behind=true;
      for(int k=0;k<vDepth;++k)if(inward[k][i])behind=false;
      cage.rigid[i]=behind;
    }
    if(cage.vActive[i]&&cage.rigid[i])
      throw std::runtime_error("Velocity active cage overlaps rigid interior");
    if(cage.vActive[i] && (bool(cage.interior[i])!=vInside))
      throw std::runtime_error("Velocity cage side mismatch");
    if(cage.pActive[i] && (bool(cage.interior[i])!=pInside))
      throw std::runtime_error("Pressure cage side mismatch");
  }

  std::vector<int> degree(N,1);
  ogsGatherScatter(degree.data(),ogsInt,ogsAdd,mesh->ogs);
  double cntLocal[9]={0},cntGlobal[9]={0};
  double minV=1e100,maxV=0,minP=1e100,maxP=0;
  for(dlong i=0;i<N;++i){
    const double w=1.0/degree[i];
    cntLocal[0]+=w*cage.interior[i];
    for(int k=0;k<3;++k){
      cntLocal[1+k]+=w*outward[k][i];
      cntLocal[4+k]+=w*inward[k][i];
    }
    cntLocal[7]+=w*cage.vActive[i];
    cntLocal[8]+=w*cage.pActive[i];
    const double dx=x[i]-part.x,dy=y[i]-part.y,dz=z[i]-part.z;
    const double rr=std::sqrt(dx*dx+dy*dy+dz*dz)/part.radius;
    if(cage.vActive[i]){minV=std::min(minV,rr);maxV=std::max(maxV,rr);}
    if(cage.pActive[i]){minP=std::min(minP,rr);maxP=std::max(maxP,rr);}
  }
  MPI_Allreduce(cntLocal,cntGlobal,9,MPI_DOUBLE,MPI_SUM,comm);
  double globalMinV=0,globalMaxV=0,globalMinP=0,globalMaxP=0;
  MPI_Allreduce(&minV,&globalMinV,1,MPI_DOUBLE,MPI_MIN,comm);
  MPI_Allreduce(&maxV,&globalMaxV,1,MPI_DOUBLE,MPI_MAX,comm);
  MPI_Allreduce(&minP,&globalMinP,1,MPI_DOUBLE,MPI_MIN,comm);
  MPI_Allreduce(&maxP,&globalMaxP,1,MPI_DOUBLE,MPI_MAX,comm);
  const long long expected[9]={
    493,294,378,470,218,150,98,
    294,218
  };
  for(int j=0;j<9;++j)
    if(std::llround(cntGlobal[j])!=expected[j])
      throw std::runtime_error("Topology sweep cage count mismatch (offline vs runtime)");
  if(!vInside && qOverA<=globalMaxV)
    throw std::runtime_error("Sampling radius not outside active Uout cage");
  if(rank==0) {
    std::printf("[CAGE100] %s / %s counts V=%lld P=%lld; V r=[%.9g,%.9g] P r=[%.9g,%.9g]\n",
      vel.c_str(),prs.c_str(),expected[7],expected[8],
      globalMinV,globalMaxV,globalMinP,globalMaxP);
    std::ofstream f("cage_stage2_summary.csv",std::ios::trunc);
    f<<"velocity_cage,pressure_cage,global_inside_nodes,"
        "Uout_S1,Uout_S2,Uout_S3,Pin_S1,Pin_S2,Pin_S3,"
        "active_velocity_nodes,active_pressure_nodes,"
        "min_velocity_r_over_a,max_velocity_r_over_a,"
        "min_pressure_r_over_a,max_pressure_r_over_a,q_over_a\n";
    f<<std::setprecision(17)<<vel<<","<<prs;
    for(int j=0;j<9;++j)f<<","<<expected[j];
    f<<","<<globalMinV<<","<<globalMaxV<<","<<globalMinP<<","<<globalMaxP<<","<<qOverA<<"\n";
  }
  stage2Cages.clear();
  stage2Cages.push_back(std::move(cage));
}
