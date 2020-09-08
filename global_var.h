#include "global_ex.h"
#ifndef PEB_STRUCT
#define PEB_STRUCT
#define MAXNUM 50
extern double dt_ring[ring_num];
extern double dtco_ring[ring_num];
extern double sigma_pre[ring_num];
extern double sigma_aft[ring_num];

typedef struct DUST_STRUCT{
  double r;
  double rf;
  double sigma;
  double a_p;
  double St;//Stokes number
  double Nd;
  double m_peb;
  double h;
  double d_v;
  double vr;
  double f_m;
  double a_frag;
  double a_drift;
  double a_df;
  double a_gr;
  double tau_gr;
  double St0;
  double vr0;
} DUST_STRUCT;

extern DUST_STRUCT dust[ring_num+1];

typedef struct INPUT_STRUCT{
  double r;
  double temp;
} INPUT_STRUCT;

extern INPUT_STRUCT input[N_in];

typedef struct DISK_STRUCT{
  double r;
  double sigma;
  double h;
  double temp;
  double rho;
  double cs;
  double visc_mol;
  double yeta;
  double yetavk;
} DISK_STRUCT;
extern DISK_STRUCT disk[ring_num];


typedef struct PEBBLE_MAP{
        double dr;
	double rad;//inner edge radius
	double rad_med;//middle radius
	double AREA;//ring AREA
  double surf_dens_gas;//gas surface density
	double time;
        double size[peb_size_num+1];
	double size_med[peb_size_num];
	double mass_in[peb_size_num];
	double mass_out[peb_size_num];
        double surf_dens[peb_size_num];
        
				double surf_dens_impose[peb_size_num];
	double rho[peb_size_num];
	double tau_fric[peb_size_num];
	double vr[peb_size_num];//vr when both r and size are in center of box
	double vt[peb_size_num];
	double vr_med_r[peb_size_num+1];//vr when r=r_med
	double vr_med_s[peb_size_num];//vr when a_pb=a_pb_med at inner edge radius
	double vt_med_r[peb_size_num+1];
	double vr_drag[peb_size_num];
	double hei[peb_size_num];
	double fluxL[peb_size_num];
	double fluxR[peb_size_num];
	double fluxRb[peb_size_num];//outward flux from previous step
} PEBBLE_MAP;

typedef struct DUST_MAP{
        double dr;
        double rad;//inner edge radius
	double rad_med;
	double hei;
	double AREA;
	double mass_in;
	double mass_out;
        double surf_dens;
	double rho;
} DUST_MAP;


typedef struct SPLINE
{

 double begin_k2;
 float x[MAXNUM+1];
 float y[MAXNUM+1];
 int point_num;

 float end_k2;

 float k1[MAXNUM+1];
 float k2[MAXNUM+1];

 float a3[MAXNUM],a1[MAXNUM];
 float b3[MAXNUM],b1[MAXNUM];

}SPLINE,*pSPLINE;

extern double alpha;
extern double mdot;
extern double opa;
extern int ITER;
extern SPLINE opa_line;
extern pSPLINE p_opa_line;
#endif
