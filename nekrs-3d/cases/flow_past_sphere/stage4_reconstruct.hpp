#pragma once
// Stage 4A: host-side exact analytical Lamb field reconstruction.
//
// Ported from supplied Nek5000 3D set.usr:
//   physalis_lamb_vel_at_point
//   physalis_lamb_pressure_at_point
// Fixed nonrotating spherical particle; general complex n=1..3 modes.
// NO solver constraint or penalty replacement in Stage 4A.
#include "stage3_lamb_math.hpp"
#include <limits>

namespace s4 {
using s3::C;
using s3::V;
using s3::Coeff;
using s3::ParticleState;

struct Value { V u; double p; };

inline s3::Basis safeBasis(const s3::Point& pt,int n,int m)
{
  const double sn=std::sqrt(std::max(0.0,1-pt.er[2]*pt.er[2]));
  if (sn>1.e-10) return s3::basis(pt,n,m);
  s3::Basis b{};
  const double sign=pt.er[2]>=0?1.:-1.;
  const double norm=std::sqrt((2*n+1)/(4*s3::PI)*s3::factorial(n-m)/s3::factorial(n+m));
  if(m==0){
    b.Y=C(norm*(sign>0 || n%2==0?1.:-1.),0);
  } else if(m==1){
    // P_n^1(x) = sqrt(1-x*x) P'_n(x) without Condon–Shortley.
    // At either pole, the surface gradient is the Cartesian x+iy
    // harmonic tangent vector with the corresponding P'_n(±1).
    const double D=0.5*n*(n+1)*((sign>0 || ((n+1)%2==0))?1.:-1.);
    b.Y=C(0,0);
    b.grad[0]=C(norm*D,0);
    b.grad[1]=C(0,norm*D);
    b.grad[2]=C(0,0);
  }
  for(int k=0;k<3;++k){
    const int j=(k+1)%3,l=(k+2)%3;
    b.tor[k]=pt.er[j]*b.grad[l]-pt.er[l]*b.grad[j];
  }
  return b;
}

// Combined poloidal radial/tangential amplitudes for one (n, Y_nm).
// The individual monomials follow physalis_lamb_vel_at_point, including
// removable n=2 radial limits; this ensures exact no-slip at r=a.
struct Amp { double radial=0,tangent=0; };
inline void addMonomial(Amp &a,int n,double roa,double c,double exponent)
{
  a.tangent+=c/double(n*(n+1));
  // q+2 is zero for n=2 singular branch; handled separately below.
  if (std::abs(exponent+2)>1e-12) a.radial+=c/(exponent+2);
}
inline Amp pAmp(int n,double roa)
{
  const double fn=n;
  const double a1=fn*(fn+3)/(2*(2*fn+3))*std::pow(roa,n+1);
  const double a2=.25*fn*(fn-2)*std::pow(roa,-n);
  const double a3=-.25*fn*fn*(2*fn+1)/(2*fn+3)*std::pow(roa,-(n+2));
  Amp out;
  addMonomial(out,n,roa,a1,n+1);
  if(n==2) out.radial-=.25*fn*std::pow(roa,-n);
  else addMonomial(out,n,roa,a2,-n);
  addMonomial(out,n,roa,a3,-n-2);
  return out;
}
inline Amp phiAmp(int n,double roa)
{
  const double fn=n;
  const double b1=fn*(fn+1)*std::pow(roa,n-1);
  const double b2=.5*fn*(fn-2)*(2*fn+1)*std::pow(roa,-n);
  const double b3=-.5*fn*fn*(2*fn-1)*std::pow(roa,-(n+2));
  Amp out;
  addMonomial(out,n,roa,b1,n-1);
  if(n==2)out.radial-=.5*fn*(2*fn+1)*std::pow(roa,-n);
  else addMonomial(out,n,roa,b2,-n);
  addMonomial(out,n,roa,b3,-n-2);
  return out;
}

inline Value evaluate(const V &xyz,const ParticleState &particle,
                      const std::vector<Coeff> &coeff,
                      double rho,double viscosity)
{
  if(coeff.size()!=s3::modes().size() || rho<=0 || viscosity<=0 || particle.a<=0)
    throw std::runtime_error("Invalid Stage4 reconstruction input");
  V d={xyz[0]-particle.x[0],xyz[1]-particle.x[1],xyz[2]-particle.x[2]};
  const double r=std::sqrt(s3::dot(d,d));
  if(r<1.e-12*particle.a)
    throw std::runtime_error("Lamb expansion not evaluated at particle center");
  const double roa=r/particle.a, aor=1/roa, nu=viscosity/rho;
  const V er=s3::scale(d,1/r);
  const double phi=std::atan2(d[1],d[0]);
  const s3::Point pt{xyz,er,0.,std::acos(std::max(-1.,std::min(1.,er[2]))),phi};
  Value val{particle.v,0.0};
  const auto ml=s3::modes();
  for(size_t j=0;j<ml.size();++j){
    const int n=ml[j].n,m=ml[j].m;
    const double fac=m==0?1.:2.;
    const auto b=safeBasis(pt,n,m);
    const double pf=fac*std::real(coeff[j].p*b.Y);
    const double ff=fac*std::real(coeff[j].phi*b.Y);
    const Amp ar=pAmp(n,roa), af=phiAmp(n,roa);
    const double Ap=(1-.5*n*(2*n-1)/double(n+1)*std::pow(aor,2*n+1))*std::pow(roa,n);
    const double Bp=-double(n*(2*n-1)*(2*n+1))/double(n+1)*std::pow(aor,n+1);
    val.p+=(viscosity*nu/(particle.a*particle.a))*(Ap*pf+Bp*ff);
    const double torCoef=std::pow(roa,-n-1)-std::pow(roa,n);
    for(int k=0;k<3;++k){
      const double pol=fac*std::real(
          (ar.tangent*coeff[j].p+af.tangent*coeff[j].phi)*b.grad[k]
          +torCoef*coeff[j].v*b.tor[k]);
      val.u[k]+=(nu/particle.a)*((ar.radial*pf+af.radial*ff)*er[k]+pol);
    }
  }
  return val;
}
inline std::vector<Coeff> manufacturedCoefficients()
{
  std::vector<Coeff> c(s3::modes().size());
  const auto ml=s3::modes();
  for(size_t i=0;i<c.size();++i){
    const int n=ml[i].n,m=ml[i].m;
    c[i].p=C(.08*n+.013*m,m>0?-.027*m:0);
    c[i].phi=C(-.047*n+.009*m,m>0?.021*m:0);
    c[i].v=C(.011*n-.019*m,m>0?-.017*m:0);
  }
  return c;
}
inline double noSlipManufacturedTest()
{
  const ParticleState p={0,{0,0,0},{0,0,0},1.};
  const auto c=manufacturedCoefficients();
  double e=0;
  const auto pts=s3::spherePoints(p,1.14);
  for(const auto& pt:pts){
    const auto v=evaluate(pt.er,p,c,1.,2.);
    for(int k=0;k<3;++k)e=std::max(e,std::abs(v.u[k]));
    if(!std::isfinite(v.p))throw std::runtime_error("Nonfinite surface pressure");
  }
  // Axial poles (safeBasis branch) should also satisfy no-slip.
  for(const V x : {V{0,0,1},V{0,0,-1}}){
    const auto v=evaluate(x,p,c,1.,2.);
    for(int k=0;k<3;++k)e=std::max(e,std::abs(v.u[k]));
  }
  return e;
}
inline double reconstructionConsistencyTest()
{
  const ParticleState p={0,{0,0,0},{0,0,0},1.};
  const auto c=manufacturedCoefficients();
  const auto pts=s3::spherePoints(p,1.14);
  const auto simplified=s3::manufacturedSurface(pts,c,p,1.14,1.,2.,0.);
  double err=0;
  for(size_t i=0;i<pts.size();++i){
    const auto v=evaluate(pts[i].x,p,c,1.,2.);
    err=std::max(err,std::abs(v.p-simplified[i].p));
    const V tangentFull=s3::plus(v.u,s3::scale(pts[i].er,-s3::dot(v.u,pts[i].er)));
    const V tangentSimple=s3::plus(simplified[i].u,s3::scale(pts[i].er,-s3::dot(simplified[i].u,pts[i].er)));
    for(int k=0;k<3;++k)err=std::max(err,std::abs(tangentFull[k]-tangentSimple[k]));
  }
  return err;
}

inline s3::V lambForce(const std::vector<Coeff>& c,
                      double mu,double rho)
{
  const auto ml=s3::modes();
  C p10=0,p11=0,phi10=0,phi11=0;
  for(size_t i=0;i<ml.size();++i) if(ml[i].n==1){
    if(ml[i].m==0){p10=c[i].p;phi10=c[i].phi;}
    else {p11=c[i].p;phi11=c[i].phi;}
  }
  const double n0=std::sqrt(3/(4*s3::PI));
  const double n1=std::sqrt(3/(8*s3::PI));
  const double factor=s3::PI*mu*(mu/rho);
  const C z=p10+6.*phi10,xy=p11+6.*phi11;
  return {-2*n1*xy.real()*factor,2*n1*xy.imag()*factor,n0*z.real()*factor};
}
} // namespace s4
