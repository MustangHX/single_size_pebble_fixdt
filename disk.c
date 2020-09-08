#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#include "global_ex.h"

double v_K(double r){
    return 29.8*100000.0*sqrt(m_star)/sqrt(r);
  //return 2e-7/r/sqrt(r)*r*LUNIT;
}

double w_K(double r){
    return v_K(r)/(r*LUNIT);
  //return 2e-7/r/sqrt(r);
}

double alpha_func(double r){
	  double rmintr = RTRAN-DRTRAN;
	  double rmaxtr = RTRAN+DRTRAN;
	  double dr_tran = DRTRAN;
	  double r_tran = RTRAN;
		double viscosity=alpha_init;
    return alpha_init;
    if (fabs(VISCOSITYRATIO-1.0)<0.1) return viscosity;
		if (r < rmintr) viscosity *=  VISCOSITYRATIO;
		if ((r >= rmintr) && (r <= rmaxtr)) {
    viscosity *= 0.5*(1.0-VISCOSITYRATIO)*sin((r-r_tran)*M_PI/(2.0*dr_tran))+0.5*(1.0+VISCOSITYRATIO);
		}

		double c0,c1,c2,c3,rlog,alplog,rtran1,rtran2,rtran0;
		rlog=log10(r);
		if(MDOT_INT==0) {
		rtran0=0.064;rtran1=0.205;rtran2=0.73;
		if(r<rtran0)	alplog=log10(0.09848);
		else if(r>=rtran0 && r < rtran1){
			c0=-6.7811175; c1=-12.02901; c2=-8.2440562; c3=-1.8594328;
			alplog=c0+c1*rlog+c2*rlog*rlog+c3*rlog*rlog*rlog;
		}
		else if(r>=rtran1 && r<rtran2){
			c0=-4.5536009; c1=-4.1222938; c2=-0.17861849;
			alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}else if(MDOT_INT==1){
		rtran0=0.014;rtran1=0.049;rtran2=0.232;
		if(r<rtran0)  alplog=log10(0.09848);
		else if(r>=rtran0 && r < rtran1){
		c0=-8.5526396; c1=-9.0739395; c2=-3.1434803; c3=-0.22451420;
		alplog=c0+c1*rlog+c2*rlog*rlog+c3*rlog*rlog*rlog;
													    }
		else if(r>=rtran1 && r<rtran2){
		c0=-6.6333076; c1=-4.3837696; c2=-0.37942259;
		alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}else if(MDOT_INT==2){
		rtran2=0.066;
		if(r<rtran2){
		c0=-11.604954; c1= -9.0560947;c2= -2.2069353;
		alplog=c0+c1*rlog+c2*rlog*rlog;
		}
		else alplog=-4.;
		}
		if (SINEALPHA==1) return viscosity;
		else return pow(10,alplog);
}

double temperature (double r){
//temporary
if (!ITER) {
        opa=func_line1(r,p_opa_line);
        //printf("rad=%e\tOPA=%e\n",r,opa);
}

if(readin>0){
  /*FILE *f_rad, *f_T;
  double rad_in[N_in]={0.},temp_in[N_in]={0.};
  double rad=0.,Tgas;
  f_rad=fopen("rad_in.txt","r");
  f_T=fopen("temperature_xh1.75e6_1.txt","r");
  int i,imax;
  i=0;
  while(rad<=dust[ring_num-1].r && i<N_in){
    fscanf(f_rad,"%lf",&rad);
    fscanf(f_T,"%lf",&Tgas);
    rad_in[i]=rad+1e-4;
    temp_in[i]=Tgas;
  //  printf("read in rad=%e\tTgas=%e\tidx=%d\n",\
      rad,Tgas,i);
    i=i+1;
    //printf("i=%d\tradmax=%e\n",i,dust[ring_num-1].r);
    }
  fclose(f_rad);
  fclose(f_T);*/

  int i,imax;
  imax=N_in;
  for(i=0;i<imax;i++){
    if(input[i].r>=r) break;//{printf("i=%\d\n",i); break;}
  }
  //printf("r=%e\trad_in[i]=%e\ti=%d\timax=%d\n",r,rad_in[i],i,imax); 
  //return temp_in[i]+(temp_in[i]-temp_in[i-1])/(rad_in[i]\
         -rad_in[i-1])*(r-rad_in[i-1]);
  return input[i].temp+(input[i].temp-input[i-1].temp)/(input[i].r\
         -input[i-1].r)*(r-input[i-1].r);
  }

alpha = alpha_func(r);
r=r*LUNIT;
double temper_active,temper_passive;
if (!ITER) {
	opa=func_line1(r/LUNIT,p_opa_line);
	//printf("rad=%e\tOPA=%e\n",r,opa);
}
 temper_active=pow(3.0,0.2)*pow(2.0,-1.4)*pow(M_PI,-0.4)\
        *pow(mu*m_p/gamma0/k_B,0.2)*pow(opa/sig_sb,0.2)\
	*pow(alpha,-0.2)*pow(G*m_star*MUNIT,0.3)\
	*pow((1-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.4)*pow(r,-0.9);

//temper_passive=temp0*pow(r/LUNIT,-3.0/7.0);
//temper_passive=temp0*pow(r/LUNIT,-0.5);//match MMSN
temper_passive=221.*pow(r/LUNIT,-0.5);
//if(r/LUNIT<8.) temper_passive=temp0*pow(8.,-0.5)*pow(r/8./LUNIT,-2);

//if ( (mdot<0e-10 && alpha>20e-4) || r/LUNIT>10.0 \
    || (temper_passive > temper_active && \
      r/LUNIT > 0.1)) 
return temper_passive;
//else return temper_active;
}

double sound_sp(double r) {//return c_s in cgs
if (!ITER) opa=func_line1(r,p_opa_line);
return sqrt(gamma0*k_B*temperature(r)/mu/m_p);
}

double height(double r) {// return scale height in cgs
if (!ITER) opa=func_line1(r,p_opa_line);
return sound_sp(r)/w_K(r);
}

double Sigma_gas (double r) {//return surface density in cgs
alpha = alpha_func(r);
double siggas;
if (!ITER) opa=func_line1(r,p_opa_line);
//printf("alpha=%e\t sig=%e\n",alpha,mdot*MUNIT/TUNIT/3.0/M_PI/(alpha*sound_sp(r)*height(r)));
 // siggas=mdot*MUNIT/TUNIT/3.0/M_PI/(alpha*sound_sp(r)*height(r));
  //printf("r=%e\t siggas=%e\t",r,siggas);
  //siggas=1700*pow(r,-1.5);
  siggas=4.1*pow(r/rc,0.)*atan(pow(r/rc,10));
  //zsj profile
  double amp, mu_peak, sig, floor, gaussian;
  amp = 29.5*100*0.02;
  mu_peak = 70.;
  sig = 10.;
  floor = 0.01;
  gaussian = amp * exp(-0.5*((r-mu_peak)/sig)*((r-mu_peak)/sig));
  if ((r <= mu_peak) || (gaussian > amp*floor)){
    siggas= gaussian;
  }else{
    siggas=amp*floor;}
  return siggas;
}

double density(double r) {
if (!ITER) opa=func_line1(r,p_opa_line);
return Sigma_gas(r)/height(r)/sqrt(2*M_PI);
}

double mean_path(double r){

if (!ITER) opa=func_line1(r,p_opa_line);
return mu*m_p/density(r)/2e-15;
}

double viscosity( double r){//http://www.ifu.ethz.ch/IE/education/AQAM/GasKinetics
	    //return 0.4991*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r)*100000.0;
if (!ITER) opa=func_line1(r,p_opa_line);
	return 0.5*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r);
}

double pressure(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
	return temperature(r)*k_B*density(r)/mu/m_p;
}

double k_P_func(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
        double dr,i_r,r1,r2;
        int i;
        i_r=ring_num*log(r/rmin)/log(rmax/rmin);
        i=floor(i_r);
        //if(!ITER) printf("k_P func, r=%f\ti_r=%f\n",r,i_r);
        dr=0.002;
        r1=rmin*exp(i*1.0/ring_num*log(rmax/rmin));
        r2=rmin*exp((i+1)*1.0/ring_num*log(rmax/rmin));
        //printf("press1=%e\tpress2=%e\tr1=%e\tr2=%e\n",pressure(r1),pressure(r2),r1,r2);
        //return -1.0*(log(pressure(r1))-log(pressure(r2)))/(log(r1)-log(r2));
        //printf("press1=%e\tpress2=%e\tr1=%e\tr2=%e\n",pressure(r),pressure(r+dr),r,dr);
        return -1.0*(log(pressure(r))-log(pressure(r+dr)))/(log(r)-log(r+dr));
}

double yeta(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
    //printf("k_P=%e\t r=%e\n",k_P_func(r),r);
    return k_P_func(r)*(sound_sp(r)/v_K(r))*(sound_sp(r)/v_K(r));
    //return 0.0022*sqrt(r);
}

double vt_gas(double r){
if (!ITER) opa=func_line1(r,p_opa_line);
        return v_K(r)*sqrt(1-yeta(r));
}


double vr_gas (double r){
alpha = alpha_func(r);
if (!ITER) opa=func_line1(r,p_opa_line);
return 3.0*alpha*sound_sp(r)*height(r)/2.0/r/LUNIT;
}
