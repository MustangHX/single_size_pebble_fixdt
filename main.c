#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char *argv[]){
  int i;
  double tsum=0.,tsum0=0., cfl=0.3,dt=dt_fix,dt_ad,outtime=0., tot_mass=0.;
  char outdustsig[256], outdustsize[256], \
    outdustvr[256], outdustst[256], outmass[256],\
    dust_frag[256],dust_drift[256], dust_df[256], dust_gr[256], tau_gr[256];
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9, *fp10;
  
  int Restarting = 0, i_arg=0, NbRestart=0;
  
  mdot=mdot_init;
  init(tsum);
  check();

  for(i_arg=0; i_arg< argc; i_arg++){
    printf("%s\n",argv[i_arg]);
    if (strchr (argv[i_arg], 's')) {
      Restarting = 1;
      i_arg++;
      NbRestart = atoi(argv[i_arg]);
    } 
    if (NbRestart < 0) printf("Incorrect restart number\n");
  }
  
  fp1=fopen("rad.txt","w");
  fp2=fopen("gassig.txt","w");
  fp3=fopen("gasrho.txt","w");
  fp4=fopen("gashei.txt","w");
  fp5=fopen("gastemp.txt","w");
  fp6=fopen("yeta.txt","w");
  fp7=fopen("yetavk.txt","w");
  fp8=fopen("cs.txt","w");
  
  if(Restarting == 1){
    Restart(NbRestart);
    tsum=NbRestart*1.0;
    outtime=NbRestart*1.0;
  }


  for(i=0;i<ring_num;i++){
    fprintf(fp1,"%e\n",dust[i].r);
    fprintf(fp2,"%e\n",disk[i].sigma);
    fprintf(fp3,"%e\n",disk[i].rho);
    fprintf(fp4,"%e\n",disk[i].h);
    fprintf(fp5,"%e\n",disk[i].temp);
    fprintf(fp6,"%e\n",disk[i].yeta);
    fprintf(fp7,"%e\n",disk[i].yetavk);
    fprintf(fp8,"%e\n",disk[i].cs);
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);


  fp5=fopen("mass_check.txt","w");
  for(i=0;i<ring_num;i++){
    tot_mass+=dust[i].sigma*2*M_PI*(dust[i+1].rf*dust[i+1].rf\
        -dust[i].rf*dust[i].rf)*LUNIT*LUNIT/m_earth;
  }
  fprintf(fp5,"%e\t%e\n",0.,tot_mass);
  int t_count=0;
  if(Restarting == 1){
    t_count=floor(num_step*log(tsum/tsum0)/log(tlim/tsum0));
    printf("t_count=%d\n",t_count);
  }
  while(tsum<tlim){
    mdot=mdot_init;
    dt=dt_fix;//tsum0*(exp((t_count+1)*1.0/num_step*log(tlim/tsum0))\
       -exp(t_count*1.0/num_step*log(tlim/tsum0)));
    
    t_count++;
    //dt_ad=upwind(cfl,dt);
    dt_ad=advec_grow(cfl,dt,1);
    if (TWO_POP>0){
    grow_two_pop(dt,tsum);
    }
    //else upwind_ND(dt);
    //if(tsum>1e4) 
    diffusion(dt);
    //if( ((int)(tsum))%((int)(outp_step))==0){
    if( tsum>=outtime){
      printf("time=%e\t dt=%e\n",tsum,dt_ad);
      printf("out put\n");
      sprintf(outdustsig,"dust_sigma%d.txt",(int)outtime);
      sprintf(outdustsize,"dust_size%d.txt",(int)outtime);
      sprintf(outdustvr,"dust_vr%d.txt",(int)outtime);
      sprintf(outdustst,"dust_st%d.txt",(int)outtime);
    
      if( TWO_POP>0){
      sprintf(dust_frag,"dust_frag%d.txt",(int)outtime);
      sprintf(dust_drift,"dust_drift%d.txt",(int)outtime);
      sprintf(dust_df,"dust_df%d.txt",(int)outtime);
      sprintf(dust_gr,"dust_gr%d.txt",(int)outtime);
      sprintf(tau_gr,"tau_gr%d.txt",(int)outtime);
      }
      
      
      
      
      fp1=fopen(outdustsig,"w");
      fp2=fopen(outdustsize,"w");
      fp3=fopen(outdustvr,"w");
      fp4=fopen(outdustst,"w");
      if( TWO_POP >0){
      fp6=fopen(dust_drift,"w");
      fp7=fopen(dust_frag,"w");
      fp8=fopen(dust_df,"w");
      fp9=fopen(dust_gr,"w");
      fp10=fopen(tau_gr,"w");
      }



      for(i=0;i<ring_num;i++){
        fprintf(fp1,"%e\n",dust[i].sigma);
        fprintf(fp2,"%e\n",dust[i].a_p);
        fprintf(fp3,"%e\n",dust[i].vr);
        fprintf(fp4,"%e\n",dust[i].St);
        if(TWO_POP>0){
        fprintf(fp6,"%e\n",dust[i].a_drift);
        fprintf(fp7,"%e\n",dust[i].a_frag);
        fprintf(fp8,"%e\n",dust[i].a_df);
        fprintf(fp9,"%e\n",dust[i].a_gr);
        fprintf(fp10,"%e\n",dust[i].tau_gr);
        }

      }
      tot_mass=0.;
      fp5=fopen("mass_check.txt","a+");
      for(i=0;i<ring_num;i++){
        tot_mass+=dust[i].sigma*2*M_PI*(dust[i+1].rf*dust[i+1].rf\
            -dust[i].rf*dust[i].rf)*LUNIT*LUNIT/m_earth;
      }
      fprintf(fp5,"%e\t%e\n",tsum,tot_mass);
      printf("%e\t%e\n",tsum,tot_mass);


      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
      if(TWO_POP > 0){
      fclose(fp6);
      fclose(fp7);
      fclose(fp8);
      fclose(fp9);
      fclose(fp10);
      }


      outtime+=outp_step;

    }
    
    tsum+=dt;

  }
  return 0;
}
