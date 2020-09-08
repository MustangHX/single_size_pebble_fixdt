#include "global_var.h"

DUST_STRUCT dust[ring_num+1];
DISK_STRUCT disk[ring_num];
INPUT_STRUCT input[N_in];
double dt_ring[ring_num];
double dtco_ring[ring_num];
double sigma_pre[ring_num];
double sigma_aft[ring_num];
double alpha;
double mdot;
double opa;
int ITER;
SPLINE opa_line;

pSPLINE p_opa_line = &opa_line;
