#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

double advec_grow(double cfl, double dt0, int GROW){
  int i;
  double dt,umax=0.,u,u01,u1,u2, uf, uf1, qx,r,r01,r1,r2,\
rf,rf1,rf01,St,St1,St2, St01,Stf,Stf1,a_p,a_p01,a_p1, a_p2; //r01 is i-1,r1 is i+1
  double sigtemp[ring_num]={0.},Ndtemp[ring_num]={0.},\
         Nd_arr[ring_num]={0.}, aptemp[ring_num]={0.};
  double d_sig,d_Nd,d_Ndgr;
  double srcterm, srcterm0, srcterm1, d_v, hdust, grow_fac, N_d;
  
  for(i=0;i<ring_num;i++){
    sigtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
    Ndtemp[i]=dust[i].r*LUNIT*dust[i].Nd;
    Nd_arr[i]=dust[i].Nd;
    aptemp[i]=dust[i].a_p;
  }

  for(i=0;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;

    a_p=aptemp[i];
    a_p01=aptemp[i-1];
    a_p1=aptemp[i+1];

    St=stokes(r/LUNIT,a_p);
    Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
    Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
    
    u=v_r(r/LUNIT,St);
    uf=v_r(rf/LUNIT,Stf);
    uf1=v_r(rf1/LUNIT,Stf1);
    

    d_Nd=(uf*Ndtemp[i]-uf1*Ndtemp[i+1])/(rf1-rf);
    
    
    if (!(GROW>0) && fabs(d_Nd)/Ndtemp[i]>umax){
      umax=fabs(d_Nd)/Ndtemp[i];
    }

    if(GROW>0){
      N_d=Nd_arr[i];
      d_v=v_pp(r/LUNIT,a_p,0);
      alpha=alpha_func(r/LUNIT);
      hdust=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
      srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
      if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
      grow_fac=1.0;
      if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
      d_Ndgr=grow_fac*srcterm;
      if(fabs(0.*d_Nd/r+d_Ndgr)>umax){
        umax=fabs(0.*d_Nd/r+d_Ndgr)/Nd_arr[i];
      }
    }

  }
  dt=cfl/umax;
  //printf("dt=%e\n",dt/TUNIT);
  dt=dt0*TUNIT;
  for(i=0;i<ring_num;i++){
    sigtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
    Ndtemp[i]=dust[i].r*LUNIT*dust[i].Nd;
  }
  
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;

    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    
    a_p=aptemp[i];
    a_p01=aptemp[i-1];
    a_p1=aptemp[i+1];

    
    //St=dust[i].St;
    St=stokes(r/LUNIT,a_p);
    //if(fabs(St1-dust[i+1].St)/St1>0.1) printf("St1_0=%e\tSt1=%e\tr1=%e\ta_p=%e\n",dust[i+1].St,St1,r1/LUNIT,a_p);
    Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
    Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
    
    u=v_r(r/LUNIT,St);
    uf=v_r(rf/LUNIT,Stf);
    uf1=v_r(rf1/LUNIT,Stf1);
    
    if(TWO_POP>0){
     uf=uf*dust[i].f_m+dust[i].vr0*(1.-dust[i].f_m);
     uf1=uf1*dust[i+1].f_m+dust[i+1].vr0*(1.-dust[i+1].f_m);
    }

    if (uf>=0. && uf1>=0.){
  //    qx=(r*u*dust[i].sigma-r1*v_r(r1/LUNIT,a_p1)\
         *dust[i-1].sigma)/(r-r1);
      d_sig=(uf*sigtemp[i-1]-uf1*sigtemp[i])/(rf1-rf);
      d_Nd=(uf*Ndtemp[i-1]-uf1*Ndtemp[i])/(rf1-rf);
      
      //printf("outward\n");
    }
    else if (uf<0. && uf1<0.){
      d_sig=(uf*sigtemp[i]-uf1*sigtemp[i+1])/(rf1-rf);
      d_Nd=(uf*Ndtemp[i]-uf1*Ndtemp[i+1])/(rf1-rf);
    //printf("qx=%e\t%e\t%e\t%e\t%e\n",qx,u,qtemp[i],rf1/LUNIT,rf/LUNIT);
    }
    else if (uf > 0.0 && uf1 < 0.0){
      d_sig=(uf*sigtemp[i-1]-uf1*sigtemp[i+1])/(rf1-rf);
      d_Nd=(uf*Ndtemp[i-1]-uf1*Ndtemp[i+1])/(rf1-rf);
    }
    else {   //(uf<0.0 && uf1 > 0.0){
      d_sig=(uf*sigtemp[i]-uf1*sigtemp[i])/(rf1-rf);
      d_Nd=(uf*Ndtemp[i]-uf1*Ndtemp[i])/(rf1-rf);
    }

    if((sigtemp[i]+dt*d_sig)/r>sigdust_floor){
      dust[i].sigma=(sigtemp[i]+dt*d_sig)/r;
    }
    else dust[i].sigma=sigdust_floor;
   
    if((Ndtemp[i]+dt*d_Nd)/r>ND_floor){
      dust[i].Nd=(Ndtemp[i]+dt*d_Nd)/r;
    }
    else dust[i].Nd=ND_floor;
    
    if(!(GROW>0)){
      dust[i].m_peb=dust[i].sigma/dust[i].Nd;
      dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
      dust[i].St=stokes(r/LUNIT,a_p);
      dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
      dust[i].d_v=d_v;
      dust[i].vr=v_r(r/LUNIT,St);
    }

    if(GROW>0){
      N_d=Nd_arr[i];
      d_v=v_pp(r/LUNIT,a_p,0);
      alpha=alpha_func(r/LUNIT);
      hdust=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
      srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
      if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
      grow_fac=1.0;
      if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
      d_Ndgr=grow_fac*srcterm;
      dust[i].Nd+=grow_fac*srcterm*dt;
      dust[i].m_peb=dust[i].sigma/dust[i].Nd;
      dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
          St=stokes(r/LUNIT,a_p);
      dust[i].St=stokes(r/LUNIT,a_p);
      dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
      dust[i].d_v=d_v;
      dust[i].vr=v_r(r/LUNIT,St);
    }

  }
  

  //boundary
  //inner 
  r=dust[0].r*LUNIT;
  rf=dust[0].rf*LUNIT;
  rf1=dust[1].rf*LUNIT;

  a_p=aptemp[0];
  a_p01=aptemp[0];
  a_p1=aptemp[1];
     
  St=stokes(r/LUNIT,a_p);
  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);


  uf=v_r(rf/LUNIT,Stf);
  uf1=v_r(rf1/LUNIT,Stf1);
  if(TWO_POP>0){
    uf=uf*dust[0].f_m+dust[0].vr0*(1.-dust[0].f_m);
    uf1=uf1*dust[1].f_m+dust[1].vr0*(1.-dust[1].f_m);
  }

  if(uf<0. && uf1<0.){
    d_sig=(1.*uf*sigtemp[0]-uf1*sigtemp[1])/(rf1-rf);
    d_Nd=(1.*uf*Ndtemp[0]-uf1*Ndtemp[1])/(rf1-rf);
  }
  else if(uf>0. && uf1>0.){
    d_sig=(uf*0.-uf1*sigtemp[0])/(rf1-rf);
    d_Nd=(uf*0.-uf1*Ndtemp[0])/(rf1-rf);
  }

  if((sigtemp[0]+dt*d_sig)/r>sigdust_floor){
    dust[0].sigma=(sigtemp[0]+dt*d_sig)/r;
  }
  else dust[0].sigma=sigdust_floor;

  if((Ndtemp[0]+dt*d_Nd)/r>ND_floor){
    dust[0].Nd=(Ndtemp[0]+dt*d_Nd)/r;
  }
  else dust[0].Nd=ND_floor;
  
  if(!(GROW>0)){
    dust[0].m_peb=dust[i].sigma/dust[i].Nd;
    dust[0].a_p=cbrt(3*dust[0].m_peb/4/M_PI/rho_peb);
    St=stokes(r/LUNIT,a_p);
    dust[0].St=stokes(r/LUNIT,a_p);
    dust[0].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[0].d_v=d_v;
    dust[0].vr=v_r(r/LUNIT,St);
  }
  if(GROW>0){
    N_d=Nd_arr[0];
    d_v=v_pp(r/LUNIT,a_p,0);
    alpha=alpha_func(r/LUNIT);
    hdust=disk[0].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
    if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
    grow_fac=1.0;
    if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    d_Ndgr=grow_fac*srcterm;
    dust[0].Nd+=d_Ndgr*dt;
    dust[0].m_peb=dust[0].sigma/dust[0].Nd;
    dust[0].a_p=cbrt(3*dust[0].m_peb/4/M_PI/rho_peb);
    St=stokes(r/LUNIT,a_p);
    dust[0].St=stokes(r/LUNIT,a_p);
    dust[0].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[0].d_v=d_v;
    dust[0].vr=v_r(r/LUNIT,St);
  }

  //outer
  r=dust[ring_num-1].r*LUNIT;
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


  if(TWO_POP>0){
    uf=uf*dust[ring_num-1].f_m+dust[ring_num-1].vr0*(1.-dust[ring_num-1].f_m);
  }
  
  d_sig=uf*(sigtemp[ring_num-1]-0.)/(rf1-rf);
  d_Nd=uf*(Ndtemp[ring_num-1]-0.)/(rf1-rf);
  
  
  if((sigtemp[ring_num-1]+dt*d_sig)/r>sigdust_floor){
    dust[ring_num-1].sigma=(sigtemp[ring_num-1]+dt*d_sig)/r;
  }
  else dust[ring_num-1].sigma=sigdust_floor;

  if((Ndtemp[ring_num-1]+dt*d_Nd)/r>ND_floor){
    dust[ring_num-1].Nd=(Ndtemp[ring_num-1]+dt*d_Nd)/r;
  }
  else dust[ring_num-1].Nd=ND_floor;
  
  if(!(GROW>0)){
    dust[ring_num-1].m_peb=dust[i].sigma/dust[i].Nd;
    dust[ring_num-1].a_p=cbrt(3*dust[ring_num-1].m_peb/4/M_PI/rho_peb);
    St=stokes(r/LUNIT,a_p);
    dust[ring_num-1].St=stokes(r/LUNIT,a_p);
    dust[ring_num-1].h=disk[ring_num-1].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[ring_num-1].d_v=d_v;
    dust[ring_num-1].vr=v_r(r/LUNIT,St);
  }
  if(GROW>0){
    N_d=Nd_arr[ring_num-1];
    d_v=v_pp(r/LUNIT,a_p,0);
    alpha=alpha_func(r/LUNIT);
    hdust=disk[ring_num-1].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));
    srcterm=-2*sqrt(M_PI)*a_p*a_p*N_d*N_d*d_v/hdust;
    if (a_p<0.) printf("src%e\ta_p%e\tN_d%e\td_v%e\thdust%e\tr%e\n",srcterm,a_p,N_d,d_v,hdust,r/LUNIT);
    grow_fac=1.0;
    if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    d_Ndgr=grow_fac*srcterm;
    dust[ring_num-1].Nd+=d_Ndgr*dt;
    dust[ring_num-1].m_peb=dust[ring_num-1].sigma/dust[ring_num-1].Nd;
    dust[ring_num-1].a_p=cbrt(3*dust[ring_num-1].m_peb/4/M_PI/rho_peb);
    St=stokes(r/LUNIT,a_p);
    dust[ring_num-1].St=stokes(r/LUNIT,a_p);
    dust[ring_num-1].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[ring_num-1].d_v=d_v;
    dust[ring_num-1].vr=v_r(r/LUNIT,St);
  }

  return dt/TUNIT;
}


double simple_advection(double cfl){
  int i;
  double dt,umax=0.,u,u1,u2,qx,r,r1,r2,St,St1,St2,Stf,Stf1,a_p,a_p01,a_p1; //r01 is i-1,r1 is i+1
  double qtemp[ring_num]={0.};
  for(i=0;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].St)/(r2-r))>umax){
      umax=fabs(v_r(dust[i].r,dust[i].St))/(r2-r);
    }
  }
  dt=cfl/umax;
  //dt=1e1*TUNIT;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i-1].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;
    St=dust[i].St;
    St1=dust[i-1].St;
    St2=dust[i+1].St;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    u=v_r(r/LUNIT,St);
    u1=v_r(r1/LUNIT,St1);
    u2=v_r(r2/LUNIT,St2);

    if (u>0.){
      
    }
  }

  return dt/TUNIT;
}



double upwind(double cfl, double dt0){
  int i;
  double dt,umax=0.,u,u01,u1,u2, uf, uf1, qx,r,r01,r1,r2,\
rf,rf1,rf01,St,St1,St2, St01,Stf,Stf1,a_p,a_p01,a_p1, a_p2; //r01 is i-1,r1 is i+1
  double qtemp[ring_num]={0.};
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
  dt=dt0*TUNIT*vr_fac;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r01=dust[i-1].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;

    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    St=dust[i].St;
    St01=dust[i-1].St;
    St1=dust[i+1].St;
    St1=stokes(r1/LUNIT,a_p);
    //if(fabs(St1-dust[i+1].St)/St1>0.1) printf("St1_0=%e\tSt1=%e\tr1=%e\ta_p=%e\n",dust[i+1].St,St1,r1/LUNIT,a_p);
    Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
    Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
    u=v_r(r/LUNIT,St);
    uf=v_r(rf/LUNIT,Stf);
    u1=v_r(r1/LUNIT,St1);
    uf1=v_r(rf1/LUNIT,Stf1);


//    r01=dust[i-1].r*LUNIT;
//    r1=dust[i+1].r*LUNIT;
//    a_p=dust[i].a_p;
//    St=dust[i].St;
//    St01=dust[i-1].St;
//    St1=dust[i+1].St;
//    a_p01=dust[i-1].a_p;
//    a_p1=dust[i+1].a_p;
//    u=v_r(r/LUNIT,St);
//    u01=v_r(r01/LUNIT,St01);
//    u1=v_r(r1/LUNIT,St1);

    if(TWO_POP>0){
     uf=uf*dust[i].f_m+dust[i].vr0*(1.-dust[i].f_m);
     //printf("i=%d\tu=%e\f_m=t%e\n",i,u,dust[i].f_m);
     //u01=u01*dust[i-1].f_m+dust[i-1].vr0*(1.-dust[i-1].f_m);
     uf1=uf1*dust[i+1].f_m+dust[i+1].vr0*(1.-dust[i+1].f_m);
    }

    if (uf>=0. && uf1>=0.){
  //    qx=(r*u*dust[i].sigma-r1*v_r(r1/LUNIT,a_p1)\
         *dust[i-1].sigma)/(r-r1);
      qx=(uf*qtemp[i-1]-uf1*qtemp[i])/(rf1-rf);
      //printf("outward\n");
    }
    else if (uf<0. && uf1<0.){
  //    qx=(r2*v_r(r2/LUNIT,a_p2)*dust[i+1].sigma\
         -r*u*dust[i].sigma)/(r2-r);
      qx=(uf*qtemp[i]-uf1*qtemp[i+1])/(rf1-rf);
      //qx=(u*qtemp[i]-u1*qtemp[i+1])/(r1-r);
//      qx=-1.*(u1*qtemp[i+1]-u*qtemp[i])/(r1-r);

    //printf("qx=%e\t%e\t%e\t%e\t%e\n",qx,u,qtemp[i],rf1/LUNIT,rf/LUNIT);

      //if(qx<0 && r/LUNIT>2. && r/LUNIT < 3. ) printf("qx r=%e\t",r/LUNIT);
    }
    else if (uf > 0.0 && uf1 < 0.0){
      //qx=(u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
      qx=(uf*qtemp[i-1]-uf1*qtemp[i+1])/(rf1-rf);
    }
    else {//(uf<0.0 && uf1 > 0.0){
      //qx=(-1.*u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
      qx=(uf*qtemp[i]-uf1*qtemp[i])/(rf1-rf);
    }
    //else if (u1<0. && u<0. && u2>0.){
     // qx=(-1.*u*qtemp[i])/(r-r1);
   // }
    //else if (u1<0. && u>0. && u2>0.){
     // qx=(u*qtemp[i])/(r-r1);
   // }

    if((qtemp[i]+dt*qx)/r>sigdust_floor){
      dust[i].sigma=(qtemp[i]+dt*qx)/r;
    }
    else dust[i].sigma=sigdust_floor;
  }
  //boundary 
  r=dust[0].r*LUNIT;
  r1=dust[1].r*LUNIT;
  rf=dust[0].rf*LUNIT;
  rf1=dust[1].rf*LUNIT;
  St=dust[0].St;
  St01=dust[1].St;
  u=v_r(r/LUNIT,St);
  u1=v_r(r1/LUNIT,St01);
  
  a_p=dust[0].a_p;
  a_p01=dust[0].a_p;
  a_p1=dust[1].a_p;

  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);


  uf=v_r(rf/LUNIT,Stf);
  uf1=v_r(rf1/LUNIT,Stf1);
  if(TWO_POP>0){
    uf=uf*dust[0].f_m+dust[0].vr0*(1.-dust[0].f_m);
    uf1=uf1*dust[1].f_m+dust[1].vr0*(1.-dust[1].f_m);
  }

  if(uf<0. && uf1<0.){          
    qx=(1.*uf*qtemp[0]-uf1*qtemp[1])/(rf1-rf);
  }
  else if(uf>0. && uf1>0.){
    qx=(uf*0.-uf1*qtemp[0])/(rf1-rf);
  }
  if((qtemp[0]+dt*qx)/r>sigdust_floor){
    dust[0].sigma=(qtemp[0]+dt*qx)/r;
  }
  else dust[0].sigma=sigdust_floor;

  r=dust[ring_num-1].r*LUNIT;
  r01=dust[ring_num-2].r*LUNIT;
  rf=dust[ring_num-1].rf*LUNIT;
  rf1=dust[ring_num].rf*LUNIT;
  St=dust[ring_num-1].St;

  a_p=dust[ring_num-1].a_p;
  a_p01=dust[ring_num-2].a_p;
  a_p1=dust[ring_num-1].a_p;

  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);

  uf=v_r(rf/LUNIT,Stf);
  u=v_r(r/LUNIT,St);


  if(TWO_POP>0){
    uf=uf*dust[ring_num-1].f_m+dust[ring_num-1].vr0*(1.-dust[ring_num-1].f_m);
  }
  qx=uf*(qtemp[ring_num-1]-0.)/(rf1-rf);
  if((qtemp[ring_num-1]+dt*qx)/r>sigdust_floor){
    dust[ring_num-1].sigma=(qtemp[ring_num-1]+dt*qx)/r;
    //printf("qx=%e\t%e\t%e\t%e\t%e\n",qx,u,qtemp[ring_num-1],rf1/LUNIT,rf/LUNIT);
  }
  else dust[ring_num-1].sigma=sigdust_floor;


  return dt/TUNIT/vr_fac;
}


double van_Leer(double cfl){
  double dt, vr_fac;
  return dt/TUNIT/vr_fac;
}
