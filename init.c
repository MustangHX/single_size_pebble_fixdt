#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

void opa_init(){// initialize opacity profile
  int i,nr=30;
  double rad[nr],opac1d[nr];
  for(i=nr-1;i>=0;i--){
  printf("seg??=%d\n",i);
  if (i>nr-2)  opa=0.1;
  else opa=opac1d[i+1];
//  rad[i]=(r_min*0.9)*exp(i*1.0/nr*log(R_OUT*1.5/(r_min*0.9)));
  rad[i]=rmin+(rmax-rmin)*i*1.0/nr;
  opac1d[i]=10.;//opa_iter(rad[i],opa);
  printf("OPA_SAMPLE=%g\t%g\n",rad[i],opac1d[i]);
  opa_line.x[i]=rad[i];
  printf("opalinex=%e\n",opa_line.x[i]);
  opa_line.y[i]=opac1d[i];
  printf("opaliney=%e\n",opa_line.y[i]);
  }
  opa_line.point_num=nr;
  opa_line.begin_k2=0.0;
  opa_line.end_k2=0.0;

  line1(p_opa_line);
}

void init(double tot_t){
  tot_t*=TUNIT;
  int i;
  for(i=0;i<ring_num+1;i++){
    //dust[i].r=rmin+(rmax-rmin)*i*1.0/ring_num;
    dust[i].rf=rmin*exp(i*1.0/ring_num*log(rmax/rmin));
  }
  for(i=0;i<ring_num;i++){
    dust[i].r=(dust[i].rf+dust[i+1].rf)/2.;
  }
  if(readin>0){
  FILE *f_rad, *f_T;
  double rad_in[N_in]={0.},temp_in[N_in]={0.};
  double rad=0.,Tgas;
  f_rad=fopen("r_in.txt","r");
  f_T=fopen("T_in.txt","r");
  int i,imax;
  i=0;
  while(rad<=dust[ring_num-1].r && i<N_in){
    fscanf(f_rad,"%lf",&rad);
    fscanf(f_T,"%lf",&Tgas);
    input[i].r=rad+1e-4;
    input[i].temp=Tgas;
  //  printf("read in rad=%e\tTgas=%e\tidx=%d\n",\
      rad,Tgas,i);
    i=i+1;
    //printf("i=%d\tradmax=%e\n",i,dust[ring_num-1].r);
    }
  fclose(f_rad);
  fclose(f_T);
  }

  for(i=0;i<ring_num;i++){
    if(input[i].r>77.5) break;
  }
  int i_tran=i;
  for(i=0;i<i_tran;i++){
  //  input[i].temp=input[i_tran].temp*\
            pow((input[i].r/input[i_tran].r),-0.5);
  }
  if (!(mdot<2e-10 && alpha>8e-4)){
    ITER=1;
    printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    opa_init();
  }

  ITER=0;
  for(i=0;i<ring_num;i++){
    //dust[i].r=(dust[i].rf+dust[i+1].rf)/2.;
    disk[i].sigma=Sigma_gas(dust[i].r);
    disk[i].h=height(dust[i].r);
    disk[i].temp=temperature(dust[i].r);
    disk[i].rho=disk[i].sigma/disk[i].h/sqrt(2*M_PI);//density(dust[i].r);
    disk[i].visc_mol=viscosity(dust[i].r);
    disk[i].cs=sound_sp(dust[i].r);
    disk[i].yeta=yeta(dust[i].r);
    disk[i].yetavk=yeta(dust[i].r)*v_K(dust[i].r);

    printf("r=%e\t sigma_gas=%e\n",dust[i].r,\
        Sigma_gas(dust[i].r));
    if(1 || ( dust[i].r>0. && dust[i].r<199.)) \
      dust[i].sigma=Sigma_gas(dust[i].r)*dust_gas;
    dust[i].a_p=a_min;//*pow(dust[i].r,-1.57);
    double tau_grow=1./(dust[i].sigma/disk[i].sigma*w_K(dust[i].r));
    dust[i].a_gr=a_min*exp(tot_t/tau_grow);
    dust[i].m_peb=4.*M_PI*dust[i].a_p*dust[i].a_p\
                  *dust[i].a_p*rho_peb/3.;
    dust[i].Nd=dust[i].sigma/dust[i].m_peb;
    double St,alpha;
    St=stokes(dust[i].r,dust[i].a_p);
    alpha=alpha_func(dust[i].r);
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].St=St;
    dust[i].vr=v_r(dust[i].r,St);
    dust[i].St0=stokes(dust[i].r,a_min);
    dust[i].vr=v_r(dust[i].r,dust[i].St0);
    dust[i].f_m=0.0;
    
  }
  
}

void Restart(int rnum){
  double tot_t=rnum*TUNIT;
  double rad, sig_d, a_p, r;
  double St, alpha, tau_grow;
  int i;

  FILE *f_sig, *f_size;

  char dsize[256],dsig[256];

  printf("RESTART=%d\n",rnum);
  sprintf(dsig,"dust_sigma%d.txt",rnum);
  sprintf(dsize,"dust_size%d.txt",rnum);
  f_sig=fopen(dsig,"r");
  f_size=fopen(dsize,"r");

  for(i=0;i<ring_num;i++){
    fscanf(f_sig,"%lf",&sig_d);
    fscanf(f_size,"%lf",&a_p);

    dust[i].sigma=sig_d;
    dust[i].a_p=a_p;
    tau_grow=1./(dust[i].sigma/disk[i].sigma*w_K(dust[i].r));
    dust[i].a_gr=a_min*exp(tot_t/tau_grow);

    dust[i].m_peb=4.*M_PI*dust[i].a_p*dust[i].a_p\
                  *dust[i].a_p*rho_peb/3.;
    dust[i].Nd=dust[i].sigma/dust[i].m_peb;
    
    St=stokes(dust[i].r,dust[i].a_p);
    alpha=alpha_func(dust[i].r);
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].St=St;
    dust[i].vr=v_r(dust[i].r,St);
    dust[i].St0=stokes(dust[i].r,a_min);
    dust[i].vr=v_r(dust[i].r,dust[i].St0);
    dust[i].f_m=0.0;
  
  }
}
