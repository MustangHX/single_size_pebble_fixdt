#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#include "global_ex.h"

void diffusion(double dt){
  dt=dt*TUNIT;
  int i;
  double alpha,alpha01,alpha1,alphaf,alphaf1, qx,r,r1,r2,rf,rf1,\
         r01,r02,a_p,a_p01,a_p1,St,St01,St1,Stf,Stf1; //r01 is i-1,r1 is i+1
  double qtemp[ring_num]={0.}, q[ring_num]={0.}, ntemp[ring_num]={0.};
  double diffco=0.01,diffco01,diffco1,Fdiff01,Fdiff1,diffcof,diffcof1;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].sigma/disk[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;
    
    
    
    r01=dust[i-1].r*LUNIT;
    //r02=dust[i-2].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    //r2=dust[i+2].r*LUNIT;

    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    St=stokes(r/LUNIT,a_p);
    St01=stokes(r01/LUNIT,a_p01);
    St1=stokes(r1/LUNIT,a_p1);
    
    
    alpha=alpha_func(r/LUNIT);
    alpha01=alpha_func(r01/LUNIT);
    alpha1=alpha_func(r1/LUNIT);


    diffco=alpha*disk[i].cs*disk[i].h/(1+St*St);
//    a_p01=dust[i-1].a_p;
    Stf=stokes(rf/LUNIT,(a_p01+a_p)/2.);//face value
    alphaf=alpha_func(rf/LUNIT);
    diffco01=alpha01*disk[i-1].cs*disk[i-1].h/(1+St01*St01);
    diffcof=alphaf*sound_sp(rf/LUNIT)*height(rf/LUNIT)/(1+Stf*Stf);
    a_p1=dust[i+1].a_p;
    Stf1=stokes(rf1/LUNIT,(a_p1+a_p)/2.);
    alphaf1=alpha_func(rf1/LUNIT);
    diffco1=alpha1*disk[i+1].cs*disk[i+1].h/(1+St1*St1);
    diffcof1=alphaf1*sound_sp(rf1/LUNIT)*height(rf1/LUNIT)/(1+Stf1*Stf1);




    Fdiff01=diffcof*Sigma_gas(rf/LUNIT)*(qtemp[i]-qtemp[i-1])\
           /(r-r01);
    Fdiff1=diffcof1*Sigma_gas(rf1/LUNIT)*(qtemp[i+1]-qtemp[i])\
           /(r1-r);
//      Fdiff01=Sigma_gas(rf/LUNIT)*(diffco*qtemp[i]-diffco01*qtemp[i-1])\
           /(r-r01);
//      Fdiff1=Sigma_gas(rf1/LUNIT)*(diffco1*qtemp[i+1]-diffco*qtemp[i])\
           /(r1-r);


    if(Fdiff01 >= 0. && Fdiff1 >=0. ){
      ntemp[i]=dust[i].Nd+dt*(Fdiff1*rf1/dust[i+1].m_peb\
               -Fdiff01*rf/dust[i].m_peb)/r/(rf1-rf);
    }else if( Fdiff01>=0. && Fdiff1<0. ){
      ntemp[i]=dust[i].Nd+dt*(Fdiff1*rf1/dust[i].m_peb\
               -Fdiff01*rf/dust[i].m_peb)/r/(rf1-rf);
    }else if( Fdiff01<0. && Fdiff1>=0. ){
      ntemp[i]=dust[i].Nd+dt*(Fdiff1*rf1/dust[i+1].m_peb\
               -Fdiff01*rf/dust[i-1].m_peb)/r/(rf1-rf);
    }else{//Fdiff01<0. && Fdiff1<0.
      ntemp[i]=dust[i].Nd+dt*(Fdiff1*rf1/dust[i].m_peb\
               -Fdiff01*rf/dust[i-1].m_peb)/r/(rf1-rf);
    }

    //printf("diff %e\tr=%e\t",(a_p)/2.,dust[i].a_p);
    //rf1=r1;
    //rf=r01;
    q[i]=qtemp[i]*disk[i].sigma+1./r*dt\
         /(rf1-rf)*(rf1*Fdiff1-rf*Fdiff01);
    //q[i]=qtemp[i]+1./r*dt*diffco*2./(r2-r1)*((qtemp[i+1]\
         -qtemp[i])/(r2-r)*r2-\
        (qtemp[i]-qtemp[i-1])/(r-r1)*r1);
    q[i]=q[i]/disk[i].sigma;
  }
//  q[1]=q[2];
  q[0]=q[1];
  ntemp[0]=q[0]*disk[0].sigma/dust[0].m_peb;
//  q[ring_num-2]=q[ring_num-3];
  q[ring_num-1]=q[ring_num-2];
  ntemp[ring_num-1]=q[ring_num-1]*\
      disk[ring_num-1].sigma/dust[ring_num-1].m_peb;
 
  for(i=0;i<ring_num;i++){
    dust[i].sigma=q[i]*disk[i].sigma;
  //  printf("r=%e\t,mp=%e\t,Nd=%e\n",dust[i].r,dust[i].m_peb,dust[i].sigma/dust[i].m_peb);
    dust[i].Nd=ntemp[i];
    //dust[i].Nd=dust[i].sigma/dust[i].m_peb;
    dust[i].m_peb=dust[i].sigma/dust[i].Nd;
    dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
  }
}


