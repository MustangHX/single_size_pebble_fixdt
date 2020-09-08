#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double stokes(double r,double a_p){
 double t_stop, tau;
 double rho_gas=density(r);
 double lambda=mean_path(r);
 double v_thermal=sqrt(8*k_B*temperature(r)/M_PI/mu/m_p);
 //printf("r=%e\t rhogas=%e\t meanpath=%e\t vthermal=%e\n",r,rho_gas,lambda,v_thermal);
 if (a_p<2.25*lambda){
  // t_stop=rho_peb*a_p/rho_gas/v_thermal;
  tau=M_PI/2.*rho_peb*a_p/Sigma_gas(r);
 }
 else{
  // t_stop = 4.0*rho_peb*a_p*a_p/9.0/rho_gas/v_thermal/lambda;
  tau=M_PI/2.*rho_peb*a_p/Sigma_gas(r)*4*a_p/9/lambda;
 }
  if (tau>1.) return 1.;
  else return tau;

 //return w_K(r)*t_stop;
}

double v_r(double r, double St){
  double vk;
  vk=v_K(r);
  //printf("St=%e\t yeta=%e\t r=%e\n",St,yeta(r),r);
  return -1.0*St*yeta(r)*vk/(1+St*St);
  //return (50.-r)*10.;

}

double v_t(double r, double St){//peb azimutal v respect to gas
  double vr;
  vr=v_r(r,St);
  return 0.5*St*vr;
}

double factor(double x){//0.5<x<1.
  //if((x-0.5)/0.5<1.) return (x-0.5)/0.5;
  //else return 1.;
  return tanh((x-0.5)/0.3);
}

double v_pp(double r, double a_p, int flag){
  double St,alpha;
  St=stokes(r,a_p);
  alpha=alpha_func(r);
  double d_vr,d_vt,d_vz,h_dust,hfactor;
  d_vr=v_r(r,St)-v_r(r,ratio_st*St);
  d_vt=v_t(r,St)-v_t(r,ratio_st*St);
  h_dust=height(r)/sqrt(1+St*(1+2*St)/alpha/(1+St));
  //if(St>0.5) hfactor=sqrt(alpha/0.5/(1+St*St));
  //else hfactor=sqrt(alpha/St/(1+St*St));
  //if(hfactor>1.) h_dust=height(r);
  //else h_dust=height(r)*hfactor;
  d_vz=h_dust*w_K(r)*(St/(1+St)-\
      ratio_st*St/(1+ratio_st*St));
  double m_peb=4.*M_PI*rho_peb*a_p*a_p*a_p/3.;
  double temper=temperature(r);
  double d_v;//=sqrt(d_vr*d_vr+d_vt*d_vt+d_vz*d_vz);
  double v_Brown=sqrt(8*2*m_peb*temper*k_B/M_PI/m_peb/m_peb);
  double v_turb,Re_turb,D_turb;
  double stfac=sqrt(3.);//factor depends on ratio_st
  D_turb=alpha*sound_sp(r)*height(r);
  Re_turb=D_turb/viscosity(r);
  //simpler version from Ueda et al. 2018
  /*if(St<=sqrt(Re_turb)) v_turb=sqrt(alpha)\
    *sound_sp(r)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else if(St>sqrt(Re_turb) && St<1.) {
    v_turb=sqrt(3.)*sqrt(alpha*St)*sound_sp(r);
  }*/
  //complex version from Okuzumi et al. 2012
  
  if(St<0.5/sqrt(Re_turb)) v_turb=sqrt(alpha)\
    *sound_sp(r)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else if(St>=0.5/sqrt(Re_turb) && St<0.5){
    v_turb=(1.-factor(St*sqrt(Re_turb)))*sqrt(alpha)*sound_sp(r)*\
           pow(Re_turb,0.25)*St*(1-ratio_st)+\
           factor(St*sqrt(Re_turb))*\
           stfac*sqrt(alpha)*sound_sp(r)*sqrt(St);
  }
  //else if(St>=1./sqrt(Re_turb) && St < 0.6)\
    v_turb=stfac*sqrt(alpha)*sound_sp(r)*sqrt(St);
   
  
//else if(St>=0.1 && St<=10.){
    //v_turb=factor(St)*sqrt(3.*alpha)*sound_sp(r)*sqrt(St)+\
           (1.-factor(St))*sqrt(alpha)*sound_sp(r)*\
           sqrt(1/(1+St)+1/(1+ratio_st*St));
 // }
  
  else v_turb=sqrt(alpha)*sound_sp(r)*sqrt(1/(1+St)+1/(1+ratio_st*St));
  //if(St>0.1) v_turb=0.;
  d_v=sqrt(d_vr*d_vr+d_vt*d_vt+d_vz*d_vz\
      +v_Brown*v_Brown+v_turb*v_turb);
  if(flag == 0) return d_v;
  else if(flag == 1) return v_Brown;
  else if(flag == 2) return fabs(d_vr);
  else if(flag == 3) return fabs(d_vt);
  else if(flag == 4) return fabs(d_vz);
  else if(flag == 5) return v_turb;
  else return d_v;
}
