#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

void upwind_ND(double dt0){
  int i;
  double dt=dt0*TUNIT,umax=0.,u,u01,u1,u2, uf, uf1, qx,r,r01,r1,r2,\
rf,rf1,rf01,St,St1,St2, St01,Stf,Stf1,a_p,a_p01,a_p1, a_p2, N_d; //r01 is i-1,r1 is i+1
  double d_v, alpha, hdust, srcterm, grow_fac=1.0;
  double qtemp[ring_num]={0.}, aptemp[ring_num]={0.};
  double vr_fac=1.0;
  for(i=0;i<ring_num-1;i++){
    //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].St)/(rf1-rf))>umax){
      //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
      umax=fabs(v_r(dust[i].r,dust[i].St))/(rf1-rf);
    }
  }
  //dt=cfl/umax;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].Nd;
    aptemp[i]=dust[i].a_p;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    //r01=dust[i-1].r*LUNIT;
    //r1=dust[i+1].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;
    
    N_d=dust[i].Nd;
    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    St=dust[i].St;
    //St01=dust[i-1].St;
    //St1=dust[i+1].St;
    //St1=stokes(r1/LUNIT,a_p);
    //St=stokes(r/LUNIT,a_p);


    //if(fabs(St1-dust[i+1].St)/St1>0.1) printf("St1_0=%e\tSt1=%e\tr1=%e\ta_p=%e\n",dust[i+1].St,St1,r1/LUNIT,a_p);
    Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
    Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
    //u=v_r(r/LUNIT,St);
    uf=v_r(rf/LUNIT,Stf);
    //u1=v_r(r1/LUNIT,St1);
    uf1=v_r(rf1/LUNIT,Stf1);

        if (uf>=0. && uf1>=0.){
      qx=(uf*qtemp[i-1]-uf1*qtemp[i])/(rf1-rf);
    }
    else if (uf<0. && uf1<0.){
      qx=(uf*qtemp[i]-uf1*qtemp[i+1])/(rf1-rf);
    }
    else if (uf > 0.0 && uf1 < 0.0){
      qx=(uf*qtemp[i-1]-uf1*qtemp[i+1])/(rf1-rf);
    }
    else {//(uf<0.0 && uf1 > 0.0){
      qx=(uf*qtemp[i]-uf1*qtemp[i])/(rf1-rf);
    }
    if((qtemp[i]+dt*qx)/r>ND_floor){
       dust[i].Nd=(qtemp[i]+dt*qx)/r;
    }
    else dust[i].Nd=ND_floor;
  //src term

    d_v=v_pp(r/LUNIT,a_p,0);
    alpha=alpha_func(r/LUNIT);
    hdust=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
    if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
    grow_fac=1.0;
    if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);

   //rk2 integrator
/*
    double srcterm0, srcterm1;
    double rk2=1.0;
    double y0=a_p,y10,y1,Nd_t,Nd_t1,d_v10,hdust10;
    d_v=v_pp(r/LUNIT,a_p,0);
    alpha=alpha_func(r/LUNIT);
    hdust=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm0=-2*sqrt(M_PI)*y0*y0*N_d*N_d*d_v/hdust;
    grow_fac=1.0;
    if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    Nd_t=N_d+grow_fac*srcterm0*dt*rk2;
    y10=cbrt(3*dust[i].sigma/Nd_t/4/M_PI/rho_peb);
    d_v10=v_pp(r/LUNIT,y10,0);
    St=stokes(r/LUNIT,y10);
    hdust10=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm1=-2*sqrt(M_PI)*y10*y10*Nd_t*Nd_t*d_v10/hdust10;
    double grow_fac1=1.0;
    if(v_frag>0. && log(v_frag/d_v10)/log(5.)<1.) grow_fac1=log(v_frag/d_v10)/log(5.);
    Nd_t1=dust[i].Nd+grow_fac*srcterm0*dt*(1.-1./2./rk2)+grow_fac1*1./2./rk2*srcterm1*dt;
   */  
  //end of rk2    

    dust[i].Nd+=grow_fac*srcterm*dt;
    dust[i].m_peb=dust[i].sigma/dust[i].Nd;
    dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
    //St=stokes(r/LUNIT,a_p);
    dust[i].St=stokes(r/LUNIT,a_p);
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].d_v=d_v;
    dust[i].vr=v_r(r/LUNIT,St);
    //printf("src%e\tmpeb%e\ta_p%e\thd%e\n",srcterm,dust[i].m_peb,\
          dust[i].a_p,hdust);

  }
  //boundary 
  N_d=dust[0].Nd; 
  r=dust[0].r*LUNIT;
  r1=dust[1].r*LUNIT;
  rf=dust[0].rf*LUNIT;
  rf1=dust[1].rf*LUNIT;
  St=dust[0].St;
  St01=dust[1].St;
  u=v_r(r/LUNIT,St);
  u1=v_r(r1/LUNIT,St01);

  a_p=aptemp[0];
  a_p01=aptemp[0];
  a_p1=aptemp[1];

  St=stokes(r/LUNIT,a_p);
  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
 

  uf=v_r(rf/LUNIT,Stf);
  uf1=v_r(rf1/LUNIT,Stf1);

  if(uf<0. && uf1<0.){
    qx=(1.*uf*qtemp[0]-uf1*qtemp[1])/(rf1-rf);
  }
  else if(uf>0. && uf1>0.){
    qx=(uf*0.-uf1*qtemp[0])/(rf1-rf);
  }
  
  if((qtemp[0]+dt*qx)/r>ND_floor){
    dust[0].Nd=(qtemp[0]+dt*qx)/r;
  } 
  else dust[0].Nd=ND_floor;
  //src term

  d_v=v_pp(r/LUNIT,a_p,0);
  alpha=alpha_func(r/LUNIT);
  hdust=disk[0].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
  srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
  if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
  grow_fac=1.0;
  if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
  dust[0].Nd+=0.*grow_fac*srcterm*dt;
  dust[0].m_peb=dust[0].sigma/dust[0].Nd;
  dust[0].a_p=cbrt(3*dust[0].m_peb/4/M_PI/rho_peb);
  St=stokes(r/LUNIT,dust[0].a_p);
  dust[0].St=St;
  dust[0].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[0].d_v=d_v;
  dust[0].vr=v_r(r/LUNIT,St);
//outer boundary
  N_d=dust[ring_num-1].Nd;
  r=dust[ring_num-1].r*LUNIT;
  r01=dust[ring_num-2].r*LUNIT;
  rf=dust[ring_num-1].rf*LUNIT;
  rf1=dust[ring_num].rf*LUNIT;

  a_p=aptemp[ring_num-1];
  a_p01=aptemp[ring_num-2];
  a_p1=aptemp[ring_num-1];
  St=stokes(r/LUNIT,a_p);


  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);

  uf=v_r(rf/LUNIT,Stf);
  u=v_r(r/LUNIT,St);


  qx=uf*(qtemp[ring_num-1]-0.)/(rf1-rf);
  if((qtemp[ring_num-1]+dt*qx)/r>ND_floor){
    dust[ring_num-1].Nd=(qtemp[ring_num-1]+dt*qx)/r;
  } 
  else dust[ring_num-1].Nd=ND_floor;
  //src term
  //fwd euler


  d_v=v_pp(r/LUNIT,a_p,0);
  alpha=alpha_func(r/LUNIT);
  hdust=disk[ring_num-1].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
  srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
  if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
  grow_fac=1.0;
  if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);

  //end of fwd euler

   //rk2 integrator
  /*
  double srcterm0, srcterm1;
  double rk2=1.0;
  double y0=a_p,y10,y1,Nd_t,Nd_t1,d_v10,hdust10;
  d_v=v_pp(r/LUNIT,a_p,0);
  alpha=alpha_func(r/LUNIT);
  hdust=disk[ring_num-1].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
  srcterm0=-2*sqrt(M_PI)*y0*y0*N_d*N_d*d_v/hdust;
  grow_fac=1.0;
  if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
  Nd_t=N_d+grow_fac*srcterm0*dt*rk2;
  y10=cbrt(3*dust[ring_num-1].sigma/Nd_t/4/M_PI/rho_peb);
  d_v10=v_pp(r/LUNIT,y10,0);
  St=stokes(r/LUNIT,y10);
  hdust10=disk[ring_num-1].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
  srcterm1=-2*sqrt(M_PI)*y10*y10*Nd_t*Nd_t*d_v10/hdust10;
  double grow_fac1=1.0;
  if(v_frag>0. && log(v_frag/d_v10)/log(5.)<1.) grow_fac1=log(v_frag/d_v10)/log(5.);
  Nd_t1=dust[ring_num-1].Nd+grow_fac*srcterm0*dt*(1.-1./2./rk2)+grow_fac1*1./2./rk2*srcterm1*dt;
   */
    //end of rk2    



  dust[ring_num-1].Nd+=grow_fac*srcterm*dt;
  dust[ring_num-1].m_peb=dust[ring_num-1].sigma/dust[ring_num-1].Nd;
  dust[ring_num-1].a_p=cbrt(3*dust[ring_num-1].m_peb/4/M_PI/rho_peb);
  St=stokes(r/LUNIT,dust[ring_num].a_p);
  dust[ring_num-1].St=St;//okes(r/LUNIT,dust[ring_num].a_p);
  dust[ring_num-1].h=disk[ring_num-1].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[ring_num-1].d_v=d_v;
  dust[ring_num-1].vr=v_r(r/LUNIT,St);
    //printf("src%e\tmpeb%e\ta_p%e\thd%e\n",srcterm,dust[i].m_peb,\


}






void upwind_size(double dt0){
  double dt=dt0*TUNIT;
  int i;
  double umax=0.,u,u1,u01,qx,r,r1,r01,a_p,a_p01,a_p1; //r1 is i+1,r01 is i-1
  double qtemp[ring_num]={0.};
  double srcfac=1.0,drtfac=1.,tau_gr;
  for(i=0;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].a_p)/(r1-r))>umax){
      umax=fabs(v_r(dust[i].r,dust[i].a_p))/(r1-r);
    }
  }
  
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].m_peb;//*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r01=dust[i-1].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;

    u=v_r(r/LUNIT,a_p);
    u01=v_r(r01/LUNIT,a_p01);
    u1=v_r(r1/LUNIT,a_p1);
    
    
    
    if (u>0.){
      qx=(u*qtemp[i]-u01*qtemp[i-1])/(r-r01);
    }
    else{
      qx=(u1*qtemp[i+1]-u*qtemp[i])/(r1-r);
    }
    dust[i].m_peb=(qtemp[i]-drtfac*dt*qx);
    //a_p=dust[i].a_p;
    //printf("a_p=%e\t m_peb=%e\n",dust[i].a_p,dust[i].m_peb);

    double St,alpha;
    //a_p=cbrt(3.*dust[i].m_peb/4./M_PI/rho_peb);
    //St=stokes(r/LUNIT,a_p);
    //dust[i].h=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));

    //St=stokes(r/LUNIT,a_p);
    alpha=alpha_func(r/LUNIT);
    // growth souce term
    double d_v;
    /*double d_vr,d_vt;
    d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
    d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

    double d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
    double v_Brown=sqrt(8*qtemp[i]*disk[i].temp*k_B/M_PI/qtemp[i]/qtemp[i]);
    double v_turb,Re_turb,D_turb;
    D_turb=alpha*disk[i].cs*disk[i].h;
    Re_turb=D_turb/disk[i].visc_mol;
    if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
    else v_turb=sqrt(3.*alpha)*disk[i].cs*sqrt(St);
    if(St>0.1) v_turb=0.;
    d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);*/
    d_v=v_pp(r/LUNIT,a_p,0);
    double rho_dust=dust[i].sigma/dust[i].h/sqrt(2.*M_PI);
     
   //tau_gr=3*dust[i].m_peb*dust[i].h/2/sqrt(M_PI)/\
            a_p/a_p/Sigma_gas(dust[i].r)/dust_gas/d_v;
    double srcterm0,srcterm1;
    //general rk2 method for srcterm
    double rk2=1.0;
    double y0=a_p,y10,y1,mpebt,mpebt1,d_v10,hdust10;
    
    //srcfac=pow(r/LUNIT/100,-0.5);
    srcterm0=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;
    //printf("rhodust=%e\t src=%e\n",dust[i].a_p,d_vr); 
    mpebt=dust[i].m_peb+srcterm0*dt*rk2;
    y10=cbrt(3*mpebt/4/M_PI/rho_peb);
    d_v10=v_pp(r/LUNIT,y10,0);
    St=stokes(r/LUNIT,y10);
    hdust10=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm1=srcfac*2.*sqrt(M_PI)*y10*y10*d_v10*dust[i].sigma/hdust10;
    mpebt1=dust[i].m_peb+srcterm0*dt*(1.-1./2./rk2)+1./2./rk2*srcterm1*dt;
    y1=cbrt(3*mpebt1/4/M_PI/rho_peb);
    double grow_fac=1.0;
    if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
    //dust[i].m_peb+=srcterm*dt;
    dust[i].m_peb=mpebt1;
 
//    tau_gr=1./(dust[i].sigma/disk[i].sigma*w_K(r/LUNIT)); 
//    a_p=dust[i].a_gr*exp(dt/tau_gr);
//    dust[i].a_gr=a_p;

//    dust[i].a_p=a_p;
//    dust[i].m_peb=4./3.*M_PI*a_p*a_p*a_p*rho_peb;
    dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
    St=stokes(r/LUNIT,dust[i].a_p);
    dust[i].St=St;
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St)); 
    //dust[i].d_v=d_v10;
    dust[i].vr=v_r(dust[i].r,St);
  }
  // boundary
  //inner
  r=dust[0].r*LUNIT;
  r1=dust[1].r*LUNIT;
  a_p=dust[0].a_p;
  u=v_r(r/LUNIT,a_p);
  qx=u*(qtemp[1]-0.*qtemp[0])/(r1-r);
  double St,alpha,d_v,d_vr,d_vt,rho_dust,srcterm,grow_fac;
  St=stokes(r/LUNIT,a_p);
  alpha=alpha_func(r/LUNIT);
  dust[0].m_peb=(qtemp[0]-dt*qx);

  /*d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
  d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

  d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
  double v_Brown=sqrt(8*qtemp[0]*disk[0].temp*k_B/M_PI/qtemp[0]/qtemp[0]);
  double v_turb,Re_turb,D_turb;
  D_turb=alpha*disk[0].cs*disk[0].h;
  Re_turb=D_turb/disk[0].visc_mol;
  if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else v_turb=sqrt(3.*alpha)*disk[0].cs*sqrt(St);
  if(St>0.1) v_turb=0.;
  d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);*/
    d_v=v_pp(r/LUNIT,a_p,0);
  rho_dust=dust[0].sigma/dust[0].h/sqrt(2.*M_PI);
  srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;

  grow_fac=1.0;
  if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
  dust[0].m_peb+=srcterm*dt;
  dust[0].a_p=cbrt(3*dust[0].m_peb/4/M_PI/rho_peb);
  dust[0].St=stokes(r/LUNIT,a_p);
  dust[0].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[0].d_v=d_v;
  dust[0].vr=v_r(dust[0].r,St);
  //
  //outer
  
  r=dust[ring_num-1].r*LUNIT;
  r1=dust[ring_num-2].r*LUNIT;
  a_p=dust[ring_num-1].a_p;
  u=v_r(r/LUNIT,a_p);
  qx=u*(0.-qtemp[ring_num-1])/(r1-r);
  St=stokes(r/LUNIT,a_p);
  alpha=alpha_func(r/LUNIT);
 /* 
  d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
  d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

  d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
  v_Brown=sqrt(8*qtemp[ring_num-1]*disk[ring_num-1].temp*k_B/M_PI/qtemp[ring_num-1]/qtemp[ring_num-1]);
  D_turb=alpha*disk[ring_num-1].cs*disk[ring_num-1].h;
  Re_turb=D_turb/disk[ring_num-1].visc_mol;
  if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else v_turb=sqrt(3.*alpha)*disk[ring_num-1].cs*sqrt(St);
  if(St>0.1) v_turb=0.;
  d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);
  */
  d_v=v_pp(r/LUNIT,a_p,0);
  rho_dust=dust[ring_num-1].sigma/dust[ring_num-1].h/sqrt(2.*M_PI);
  srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;
  //dust[ring_num-1].m_peb=(qtemp[ring_num-1]-dt*qx);

  grow_fac=1.0;
  if(v_frag>0. && log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
  dust[ring_num-1].m_peb+=srcterm*dt;
  dust[ring_num-1].a_p=cbrt(3*dust[ring_num-1].m_peb/4/M_PI/rho_peb);
  St=stokes(r/LUNIT,a_p);
  dust[ring_num-1].St=St;
  dust[ring_num-1].h=disk[ring_num-1].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[ring_num-1].d_v=d_v;

}
void grow_two_pop(double dt0, double tot_time){
  double dt=dt0*TUNIT;
  double tot_t=tot_time*TUNIT;
  double tau_grow, a_frag, a_drift, a_df, v_f,\
    cs, alpha,f_m, St, St_df, a0=a_min, a1, a_p, r;
  double srcterm, srcfac=1.0, rho_dust, d_v;
  int i;
  for (i=0;i<ring_num;i++){
    r=dust[i].r;
    alpha=alpha_func(r);
    cs=disk[i].cs;
    v_f=v_frag;
    tau_grow=1./(dust[i].sigma/disk[i].sigma*w_K(r));
    a_frag=FF*(2./3./M_PI)*(disk[i].sigma/rho_peb/alpha)*(v_f*v_f/cs/cs);
    a_drift=FD*(2*dust[i].sigma/M_PI/rho_peb)*\
            (v_K(r)*v_K(r)/cs/cs)/fabs(k_P_func(r));
    St_df=v_f*v_K(r)/k_P_func(r)/cs/cs/(1.-ratio_st);
    a_df=FF*St_df*2.*disk[i].sigma/M_PI/rho_peb;

    if(a_drift < a_frag && a_drift < a_df){
      f_m= FMD;
      a1 = a_drift;
    }
    else if (a_frag < a_drift && a_frag < a_df){
      f_m = FMF;
      a1  = a_frag;
    }
    else{
      f_m = FMF;
      a1  = a_df;
    }
    /*
    a_p=dust[i].a_p;
    d_v=v_pp(r,a_p,0);
    rho_dust=dust[0].sigma/dust[0].h/sqrt(2.*M_PI);
    srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;
    dust[i].m_peb+=srcterm*dt;
    a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
*/
    //a_p=a0*exp(tot_t/tau_grow);
    a_p=dust[i].a_gr*exp(dt/tau_grow);
    dust[i].a_gr=a_p;

    if (a_p>a1) a_p=a1;
    dust[i].a_p=a_p;
    dust[i].f_m=f_m;
    dust[i].a_frag=a_frag;
    dust[i].a_drift=a_drift;
    dust[i].a_df=a_df;
    dust[i].tau_gr=tau_grow;
    St=stokes(dust[i].r,a_p);
    dust[i].St=St;
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].vr=v_r(dust[i].r,dust[i].St);
    
  }

}
