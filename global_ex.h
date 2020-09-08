#include <math.h>
#define sigdust_floor 6e-25 //dust surface density floor
#define ND_floor 1e-10 //dust surface number density floor
#define v_frag 1000. //in cm/s
#define a_min 1e-5 // in cm
#define ratio_st 0.5
#define rho_peb 1.4
//#define alpha 0.001
#define gamma0 1.0
#define m_star 1.25
//#define M_PI 3.14159265358979323846264338327950288
#define AU_km 1.49597871e8
#define LUNIT 1.49597871e13
#define MSUN 1.9891e33
#define MUNIT 1.9891e33
#define m_earth 5.97219e27
#define mdot_init 1e-9
#define alpha_init 1e-3
#define temp0 170.  //172.043610986 //temp at 1 AU
//#define r_tran 300.365  //transition from active to passive
#define TUNIT 3.15569e7
#define k_P 2.55
#define k_B 1.38e-16
#define sig_sb 5.6704e-5
//#define opa 10.0
#define r_star 0.01395141 //T_Tauri with 3 solar radius
#define mu 2.34 
#define m_p 1.66053906660e-24 //amu //1.660538921e-24 m_p
#define G 6.674e-8
#define peb_low_lim 1e-40
#define COAG_SW 0 // coagulation must be turned on to have fragmentation
#define FRAG_SW 0
#define DIFF_SW 1 //diffusion

#define TWO_POP -1
#define FF 0.37  //a_p=FF*amax in frag limited
#define FD 0.55  //a_p=FD*amax in drift limited
#define FMD 0.97 //drift limited
#define FMF 0.75 //frag limited

#define rmax 300.0
#define rmin 30.0
#define rc   65.
#define ring_num 512
#define num_step 4e5
#define dt_fix 10.0
#define tlim 2001000. //in yr
#define outp_step 1000.
#define peb_num 10
#define v_coag_max 500.0 //cm/s
#define v_tran_width 500.0 // cm/s
#define v_frag_min0  1500.0 //cm/s
#define v_frag_min1  2000.0 //cm/s
#define v_frag_max1  2000.0 //cm/s
#define v_frag_max0  2500.0 //cm/s
#define readin 1
#define N_in 512

#define frag_slope -0.2
#define interaction_min_radius 0.2 // in AU, min r with peb-peb interaction
#define size_ring 0.25
#define size_min 0.01
#define size_min_inj 1e-2
#define size_min_init 1e-2
#define outp_time 1
#define NUM_LIM 100
#define peb_size_num 90
#define peb_size_lim 0.3 //in cm
#define size_step 0.05
#define size_slope 5.0
#define size_slope_pow 1.0
#define size_slope_birn 0.0
#define dust_gas 0.01
#define peb_dust 0.01
#define peb_dust_inj 0.0
#define SINEALPHA 1
#define MDOT_INT 1 //0 for mdot=1e-8, 1 for mdot=1e-9
#define VISCOSITYRATIO 1.
#define RTRAN 50.0 
#define DRTRAN 10.0
#define PEB_IMPOSE 0
