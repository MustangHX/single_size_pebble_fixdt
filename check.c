#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

void check(){
  // v_pp at 100 au
  double a_p,r=100,r_in=10.,vr, d_v,v_Brown,\
         d_vr,d_vt,d_vz,v_turb,tgrowth,St,\
         hdust,alpha,m_peb,srcterm,N_d;
  int i,Nsize=300;
  FILE *fp,*fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8, *fp9;
  
  fp=fopen("size_check.txt","w");
  fp0=fopen("vpp_check.txt","w");
  fp1=fopen("vBrown_check.txt","w");
  fp2=fopen("dvr_check.txt","w");
  fp3=fopen("dvt_check.txt","w");
  fp4=fopen("dvz_check.txt","w");
  fp5=fopen("vturb_check.txt","w");
  fp6=fopen("tgrowth100au.txt","w");
  fp7=fopen("tgrowth_inner.txt","w");
  fp8=fopen("vr_inner.txt","w");
  fp9=fopen("srcterm_inner.txt","w");




  for(i=0;i<Nsize;i++){
    a_p=1e-5*exp(i*1.0/Nsize*log(1e3/1e-5));
    d_v=v_pp(r,a_p,0);
    v_Brown=v_pp(r,a_p,1);
    d_vr=v_pp(r,a_p,2);
    d_vt=v_pp(r,a_p,3);
    d_vz=v_pp(r,a_p,4);
    v_turb=v_pp(r,a_p,5);
    m_peb=4.*M_PI*rho_peb*a_p*a_p*a_p/3.;
    St=stokes(100.,a_p);
    alpha=alpha_func(100.);
    hdust=height(100.)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    tgrowth=3*m_peb*hdust/2/sqrt(M_PI)/\
            a_p/a_p/Sigma_gas(100.)/dust_gas/d_v;

    fprintf(fp,"%e\n",a_p);
    fprintf(fp0,"%e\n",d_v);
    fprintf(fp1,"%e\n",v_Brown);
    fprintf(fp2,"%e\n",d_vr);
    fprintf(fp3,"%e\n",d_vt);
    fprintf(fp4,"%e\n",d_vz);
    fprintf(fp5,"%e\n",v_turb);
    fprintf(fp6,"%e\n",tgrowth);
    
    d_v=v_pp(r_in,a_p,0);
    St=stokes(r_in,a_p);
    alpha=alpha_func(r_in);
    hdust=height(r_in)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    tgrowth=3*m_peb*hdust/2/sqrt(M_PI)/\
            a_p/a_p/Sigma_gas(r_in)/dust_gas/d_v;
    srcterm=2.*sqrt(M_PI)*a_p*a_p*d_v*\
            Sigma_gas(r_in)*dust_gas/hdust;
    vr=v_r(r_in,St);
    fprintf(fp7,"%e\n",tgrowth);
    fprintf(fp8,"%e\n",vr);
    fprintf(fp9,"%e\n",srcterm);

  }
  fclose(fp);
  fclose(fp0);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);
  fclose(fp9);

  double dt=1e2*TUNIT,tsum=0.0, d2g=1e-2,a_pmin=1e-5,t_int=1e4;
  fp1=fopen("fwd_eu_100yr.txt","w");
  fp2=fopen("rk2_100yr.txt","w");
  a_p=a_pmin;
  m_peb=4./3.*M_PI*a_p*a_p*a_p*rho_peb;
  while (tsum<=t_int*TUNIT){
    d_v=v_pp(r_in,a_p,0);
    St=stokes(r_in,a_p);
    alpha=alpha_func(r_in);
    hdust=height(r_in)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    /*srcterm=2.*sqrt(M_PI)*a_p*a_p*d_v*\
            Sigma_gas(r_in)*d2g/hdust;
    m_peb+=srcterm*dt;*/
    N_d=Sigma_gas(r_in)/m_peb;
    srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
    N_d+=srcterm*dt;
    m_peb=Sigma_gas(r_in)/N_d;
    a_p=cbrt(3*m_peb/4/M_PI/rho_peb);
    if(N_d<0.) printf("t=%e\tN_d=%e\tm_peb=%e\ta_p=%e\n",tsum,N_d,m_peb,a_p);
    tsum+=dt;
    fprintf(fp1,"%e\n",N_d);
  }
  fclose(fp1); 
  
/*  tsum=0.0;
  a_p=a_pmin;
  m_peb=4./3.*M_PI*a_p*a_p*a_p*rho_peb;
  double srcterm0,srcterm1;
  //general rk2 method for srcterm
  double rk2=1.0;
  double y0=a_p,y10,y1,mpebt,mpebt1,d_v10,hdust10;

  while (tsum<=t_int*TUNIT){
    d_v=v_pp(r_in,a_p,0);
    St=stokes(r_in,a_p);
    alpha=alpha_func(r_in);
    hdust=height(r_in)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    srcterm0=2.*sqrt(M_PI)*a_p*a_p*d_v*\
             Sigma_gas(r_in)*d2g/hdust;
    //printf("rhodust=%e\t src=%e\n",dust[i].a_p,d_vr); 
    mpebt=m_peb+srcterm0*dt*rk2;
    y10=cbrt(3*mpebt/4/M_PI/rho_peb);
    d_v10=v_pp(r_in,y10,0);
    St=stokes(r_in,y10);
    hdust10=height(r_in)/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm1=2.*sqrt(M_PI)*y10*y10*d_v10*\
             Sigma_gas(r_in)*d2g/hdust10;
    mpebt1=m_peb+srcterm0*dt*(1.-1./2./rk2)\
           +1./2./rk2*srcterm1*dt;
    y1=cbrt(3*mpebt1/4/M_PI/rho_peb);
    m_peb=mpebt1;
    
    a_p=cbrt(3*m_peb/4/M_PI/rho_peb);
    tsum+=dt;
    fprintf(fp2,"%e\n",a_p);
  }
  fclose(fp2);*/
}
