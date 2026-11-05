#pragma once
// Stage 3: exact m=0..n complex spherical-harmonic scalar products
// and mode-by-mode 2x2 / 1x1 coefficient extraction, based on the supplied
// NekPhysalis 3D set.usr (s36b_glphi_scalar_products, physalis_lamb_solve).
// Host-side, standalone C++17: deliberately independent of the NekRS API.
// q/a > 1; fixed, nonrotating particle; pressure frame shift = 0.
#include <array>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace s3 {
constexpr int LMAX=3, NPHI=32, NMU=32;
constexpr double PI=3.141592653589793238462643383279502884;
using C=std::complex<double>;
using V=std::array<double,3>;
using CV=std::array<C,3>;
struct ParticleState { int id; V x, v; double a; };
struct Point { V x,er; double weight,theta,phi; };
struct Fields { double p=0; V u={0,0,0}; };
struct Mode { int n,m; };
struct Coeff { C p={0,0},phi={0,0},v={0,0}; };
struct RawProducts { C sp={0,0},su={0,0},sv={0,0}; };
struct Radial {double Ap,Bp,Au,Bu,Cv;};
inline V plus(V a,V b) { return {a[0]+b[0],a[1]+b[1],a[2]+b[2]}; }
inline V scale(V a,double c) {return {c*a[0],c*a[1],c*a[2]};}
inline V cross(V a,V b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double dot(V a,V b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline double factorial(int n){double f=1; for(int k=2;k<=n;++k)f*=k;return f;}
inline std::vector<Mode> modes() {
  std::vector<Mode> out;
  for(int n=1;n<=LMAX;++n)for(int m=0;m<=n;++m)out.push_back({n,m});
  return out;
}
inline void gaussLegendre(int n,std::vector<double> &x,std::vector<double> &w) {
  x.assign(n,0);w.assign(n,0);
  for(int i=0;i<(n+1)/2;++i){
    double z=std::cos(PI*(i+0.75)/(n+0.5)), pp=0;
    bool converged=false;
    for(int j=0;j<80;++j){
      double p1=1,p2=0;
      for(int k=1;k<=n;++k){
        double p3=p2; p2=p1; p1=((2*k-1)*z*p2-(k-1)*p3)/k;
      }
      pp=n*(z*p1-p2)/(z*z-1);
      const double dz=p1/pp;
      z-=dz;
      if(std::abs(dz)<2e-15){converged=true;break;}
    }
    if(!converged)throw std::runtime_error("GL root iteration failed");
    // Re-evaluate derivative at the converged z.
    double p1=1,p2=0;
    for(int k=1;k<=n;++k) {
      const double p3=p2; p2=p1;p1=((2*k-1)*z*p2-(k-1)*p3)/k;
    }
    pp=n*(z*p1-p2)/(z*z-1);
    x[i]=-z;x[n-1-i]=z;
    const double ww=2/((1-z*z)*pp*pp);
    w[i]=ww;w[n-1-i]=ww;
  }
}
inline std::vector<Point> spherePoints(const ParticleState &p,double qOverA){
  if(!(p.a>0 && qOverA>1))throw std::runtime_error("Invalid sampling radius/particle radius");
  std::vector<double>mu,w;
  gaussLegendre(NMU,mu,w);
  std::vector<Point> pts;
  pts.reserve(NMU*NPHI);
  for(int i=0;i<NMU;++i){
    const double st=std::sqrt(std::max(0.,1-mu[i]*mu[i]));
    for(int j=0;j<NPHI;++j){
      const double phi=2*PI*j/NPHI;
      const V er={st*std::cos(phi),st*std::sin(phi),mu[i]};
      pts.push_back({plus(p.x,scale(er,p.a*qOverA)),er,
                     w[i]*2*PI/NPHI,std::acos(mu[i]),phi});
    }
  }
  return pts;
}
// No Condon–Shortley phase; matches original 3D NekPhysalis set.usr.
inline double assocP(int n,int m,double mu){
  const double st=std::sqrt(std::max(0.,1-mu*mu));
  double pmm=1;
  for(int k=1;k<=m;++k)pmm*=(2*k-1)*st;
  if(n==m)return pmm;
  double pmmp1=(2*m+1)*mu*pmm;
  if(n==m+1)return pmmp1;
  for(int k=m+2;k<=n;++k){
    const double next=((2*k-1)*mu*pmmp1-(k+m-1)*pmm)/(k-m);
    pmm=pmmp1;pmmp1=next;
  }
  return pmmp1;
}
struct Basis{C Y;CV grad, tor;};
inline Basis basis(const Point &pt,int n,int m) {
  const double mu=pt.er[2], st=std::sqrt(std::max(0.,1-mu*mu));
  if(st<1e-12)throw std::runtime_error("Stage3 GL nodes unexpectedly hit a pole");
  const double cp=std::cos(pt.phi),sp=std::sin(pt.phi);
  const V eth={mu*cp,mu*sp,-st}, eph={-sp,cp,0};
  const double pnm=assocP(n,m,mu);
  const double lower=(m<n)?assocP(n-1,m,mu):0.0;
  const double dpdth=(n*mu*pnm-(n+m)*lower)/st;
  const double norm=std::sqrt((2*n+1)/(4*PI)*factorial(n-m)/factorial(n+m));
  const C phase=std::polar(1.,double(m)*pt.phi);
  const C Y=norm*pnm*phase, dth=norm*dpdth*phase, dph=C(0,m)*Y;
  Basis b;b.Y=Y;
  for(int k=0;k<3;++k)
    b.grad[k]=dth*eth[k]+dph*(eph[k]/st);
  for(int k=0;k<3;++k){
    const int j=(k+1)%3,l=(k+2)%3;
    b.tor[k]=pt.er[j]*b.grad[l]-pt.er[l]*b.grad[j];
  }
  return b;
}
inline Radial radial(int n,double lambda){
  if(lambda<=1) throw std::runtime_error("q/a must be >1 for Lamb inversion");
  const double fn=n, roa=lambda,aor=1/lambda, inv=std::pow(aor,2*n+1);
  Radial r;
  r.Ap=(1-0.5*fn*(2*fn-1)/(fn+1)*inv)*std::pow(roa,n);
  r.Bp=-(fn*(2*fn-1)*(2*fn+1)/(fn+1))*std::pow(aor,n+1);
  r.Au=0.25*fn*(2*(fn+3)/(2*fn+3)
       +(fn-2-fn*(2*fn+1)/(2*fn+3)*aor*aor)*inv)*std::pow(roa,n+1);
  r.Bu=fn*(fn+1+0.5*((fn-2)*(2*fn+1)*roa*roa-fn*(2*fn-1))*inv)
       *std::pow(roa,n-1);
  r.Cv=fn*(fn+1)*(inv-1)*std::pow(roa,n);
  return r;
}
inline std::vector<RawProducts> project(const std::vector<Point>& pts,
                     const std::vector<Fields>& sampled,
                     const ParticleState& p,double *pMeanOut=nullptr){
  if(pts.size()!=NMU*NPHI || pts.size()!=sampled.size())
    throw std::runtime_error("Unexpected GL48x32 sample size");
  double sumP=0,sumW=0;
  for(size_t k=0;k<pts.size();++k){sumP+=pts[k].weight*sampled[k].p;sumW+=pts[k].weight;}
  const double pMean=sumP/sumW;
  if(pMeanOut)*pMeanOut=pMean;
  const auto ml=modes();
  std::vector<RawProducts> rp(ml.size());
  for(size_t k=0;k<pts.size();++k){
    V u;for(int j=0;j<3;++j)u[j]=sampled[k].u[j]-p.v[j];
    for(size_t mi=0;mi<ml.size();++mi){
      const Basis b=basis(pts[k],ml[mi].n,ml[mi].m);
      const double w=pts[k].weight;
      rp[mi].sp+=w*std::conj(b.Y)*(sampled[k].p-pMean);
      for(int d=0;d<3;++d){
        rp[mi].su+=w*std::conj(b.grad[d])*u[d];
        rp[mi].sv+=w*std::conj(b.tor[d])*u[d];
      }
    }
  }
  return rp;
}
inline std::vector<Coeff> invert(const std::vector<RawProducts>& prod,
                               const ParticleState& part,double lambda,
                               double rho,double viscosity){
  if(prod.size()!=modes().size() || rho<=0 || viscosity<=0)
    throw std::runtime_error("Invalid inversion inputs");
  const double a=part.a,nu=viscosity/rho;
  std::vector<Coeff> out(prod.size());
  const auto ml=modes();
  for(size_t i=0;i<prod.size();++i){
    const Radial r=radial(ml[i].n,lambda);
    const double det=r.Ap*r.Bu-r.Au*r.Bp;
    if(std::abs(det)<1e-14 || std::abs(r.Cv)<1e-14)
      throw std::runtime_error("Ill conditioned radial extraction system");
    const C Sp=a*a/(viscosity*nu)*prod[i].sp;
    const C Su=a/nu*prod[i].su;
    const C Sv=a/nu*prod[i].sv;
    out[i].p=(Sp*r.Bu-Su*r.Bp)/det;
    out[i].phi=(r.Ap*Su-r.Au*Sp)/det;
    out[i].v=Sv/r.Cv;
    if(ml[i].m==0){
      out[i].p=C(out[i].p.real(),0);
      out[i].phi=C(out[i].phi.real(),0);
      out[i].v=C(out[i].v.real(),0);
    }
  }
  return out;
}
// Manufactured scalar-product reconstruction at q/a; tests all complex
// p,phi,v modes. Not the full near-field radial Lamb velocity solution.
inline std::vector<Fields> manufacturedSurface(const std::vector<Point> &pts,
      const std::vector<Coeff> &coeff,const ParticleState &p,double lambda,
      double rho,double viscosity,double gauge=7.){
  if(coeff.size()!=modes().size())throw std::runtime_error("Invalid manufactured modes");
  const double nu=viscosity/rho;
  const auto ml=modes();
  std::vector<Fields> out(pts.size());
  for(size_t k=0;k<pts.size();++k){
    out[k].p=gauge;out[k].u=p.v;
    for(size_t mi=0;mi<ml.size();++mi){
      const int n=ml[mi].n,m=ml[mi].m;
      const auto r=radial(n,lambda);
      const auto b=basis(pts[k],n,m);
      const double scaleM=m>0?2.:1.;
      const C pAmp=(viscosity*nu/(p.a*p.a))*(r.Ap*coeff[mi].p+r.Bp*coeff[mi].phi);
      const C uAmp=(nu/p.a)*(r.Au*coeff[mi].p+r.Bu*coeff[mi].phi)/double(n*(n+1));
      const C vAmp=(nu/p.a)*(r.Cv*coeff[mi].v)/double(n*(n+1));
      out[k].p+=scaleM*std::real(pAmp*b.Y);
      for(int d=0;d<3;++d)
        out[k].u[d]+=scaleM*std::real(uAmp*b.grad[d]+vAmp*b.tor[d]);
    }
  }
  return out;
}
inline double manufacturedTest(){
  const ParticleState p={0,{0,0,0},{0,0,0},1};
  const auto pts=spherePoints(p,1.14);
  const auto ml=modes();
  std::vector<Coeff> expect(ml.size());
  for(size_t i=0;i<ml.size();++i){
    const int n=ml[i].n,m=ml[i].m;
    expect[i].p=C(0.08*n+0.013*m,m>0?-0.027*m:0);
    expect[i].phi=C(-0.047*n+0.009*m,m>0?0.021*m:0);
    expect[i].v=C(0.011*n-0.019*m,m>0?-0.017*m:0);
  }
  const auto synthetic=manufacturedSurface(pts,expect,p,1.14,1.,2.,7.);
  const auto prod=project(pts,synthetic,p);
  const auto actual=invert(prod,p,1.14,1.,2.);
  double err=0;
  for(size_t i=0;i<ml.size();++i){
    err=std::max({err,std::abs(expect[i].p-actual[i].p),
                     std::abs(expect[i].phi-actual[i].phi),
                     std::abs(expect[i].v-actual[i].v)});
  }
  return err;
}
}
