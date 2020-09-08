#!bin/python

import sys
sys.path.append("/Users/xiaohu/work/py_package")
from readcol import *
import matplotlib
matplotlib.use('MACOSX')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
LUNIT=1.49597871e13
TUNIT=3.15569e7
m_earth=5.97219e27
rho_peb=1.4

TWO_POP=False

name=sys.argv[1]
num=sys.argv[2]

fig=plt.figure(figsize=(10,6))
gs=gridspec.GridSpec(2,2)#,width_ratios=[1,1])

r=readcol(name+"/rad.txt")

sig0=readcol(name+"/dust_sigma0.txt")
sig=readcol(name+"/dust_sigma"+num+".txt")

size0=readcol(name+"/dust_size0.txt")
size=readcol(name+"/dust_size"+num+".txt")

if TWO_POP:
  a_drift=readcol(name+"/dust_drift"+num+".txt")
  a_frag=readcol(name+"/dust_frag"+num+".txt")
  a_df=readcol(name+"/dust_df"+num+".txt")
  a_gr=readcol(name+"/dust_gr"+num+".txt")
  tau_gr=readcol(name+"/tau_gr"+num+".txt")

rB=np.load('rB.npy')
aB=np.load('aB.npy')
sigma_dB=np.load('sigma_dB.npy')
taugrowB=np.load('taugrowB.npy')


vr0B=np.load('vr0B.npy')
vrdB=np.load('vrdB.npy')
vr0B0=np.load('vr0B0.npy')
St0B0=np.load('St0B0.npy')




vr0=readcol(name+"/dust_vr0.txt")
vr=readcol(name+"/dust_vr"+num+".txt")

St0=readcol(name+"/dust_st0.txt")
St=readcol(name+"/dust_st"+num+".txt")

ax1=plt.subplot(gs[0,0])
ax1.plot(r,sig0)
ax1.plot(r,sig)

temp_in=readcol(name+"/T_in.txt")
r_in=readcol(name+"/r_in.txt")

ax1.annotate('time='+num+" yr", xy=(0., 1.02), xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom',fontsize='large',backgroundcolor='grey',color='white')
ax1b=ax1.twinx()
rc=65.
ax1b.plot(r_in,temp_in,ls='--')#*4.1*(r/rc)**0*np.arctan((r/rc)**10),ls='--')
if TWO_POP:
  ax1.plot(rB,sigma_dB,ls='-.')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Sig_dust')
ax1.set_ylim(bottom=1e-5)
ax1.set_xlabel('r (au)')
ax1.set_xlim(30,300)
ax1b.set_ylim(0,100) 
ax1b.set_ylabel("T (K)")
ax2=plt.subplot(gs[0,1])
ax2.plot(r,size0)
ax2.plot(r,size)
if TWO_POP:
 #ax2.plot(r,a_drift,label="drift")
 #ax2.plot(r,a_frag,label="frag")
 #ax2.plot(r,a_df,label="df")
 ax2.plot(r,a_gr,label="gr")
 ax2.plot(rB,aB,ls='-.')
 ax2.set_ylim(1e-5,1e4)
 ax2.legend()


#ax2.plot(r,20*r**-1.57,ls='--',lw=1)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel('a_p (cm)')
ax2.set_xlabel('r (au)')
#ax2.set_xlim(4,100)
  
ax3=plt.subplot(gs[1,0])
ax3.plot(r,vr0)
ax3.plot(r,vr)

if TWO_POP:

  ax3.plot(rB,vr0B,ls=':')
  ax3.plot(rB,vrdB,ls='--')


ax3.set_ylabel('vr')
ax3.set_xscale('log')
#ax3.set_xlim(4,100)
ax3.set_xlabel('r (au)')

ax4=plt.subplot(gs[1,1])
#ax4.plot(r,-vr0*sig0*2*np.pi*r*LUNIT*TUNIT*1e6/m_earth)
#ax4.plot(r,-vr*sig*2*np.pi*r*LUNIT*TUNIT*1e6/m_earth)
ax4.plot(r,St0)
ax4.plot(r,St)
#ax4.plot(r,tau_gr)
#ax4.plot(rB,taugrowB,ls='-.')
#ax4.plot(rB,St0B0,ls='-.')
ax4.set_ylabel('St')
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_xlabel('r (au)')
#ax4.set_xlim(4,100)
gassig=readcol(name+"/gassig.txt")


fig.subplots_adjust(hspace=0.1,wspace=0.35,left=0.07,right=0.9,bottom=0.1, top =0.95)

#ax5=plt.subplot(gs[2,:])
#ax5.plot(r,sig0/gassig)
#ax5.plot(r,sig/gassig)
#ax5.set_yscale('log')
#ax5.set_xscale('log')
plt.show()

gastemp=readcol(name+"/gastemp.txt")

plt.plot(r,221*r**-0.5)
plt.plot(r,gastemp)
plt.xlabel("r (au)")
plt.ylabel("T (k)")
plt.show()

fig=plt.figure(figsize=(7,4))
gs=gridspec.GridSpec(1,2)

ugr=readcol(name+"/dust_ugr"+num+".txt")
uad=readcol(name+"/dust_uad"+num+".txt")

dt=1095037661.76
ax1=plt.subplot(gs[0,0])
ax1.plot(r,ugr*dt)
ax1.set_xscale('log')

ax2=plt.subplot(gs[0,1])
ax2.plot(r,uad*dt)
ax2.set_xscale('log')
plt.show()


size_check=readcol(name+"/size_check.txt")
vpp_check=readcol(name+"/vpp_check.txt");
vBrown_check=readcol(name+"/vBrown_check.txt");
dvr_check=readcol(name+"/dvr_check.txt");
dvt_check=readcol(name+"/dvt_check.txt");
dvz_check=readcol(name+"/dvz_check.txt");
vturb_check=readcol(name+"/vturb_check.txt");


fig,ax=plt.subplots()

ax.plot(size_check,vpp_check,label='vpp',c='k')
ax.plot(size_check,vBrown_check,label='vBrown',ls='--',c='grey')
ax.plot(size_check,dvr_check,label='dvr',ls='--',c='g',alpha=0.5)
ax.plot(size_check,dvt_check,label='dvt',ls='--',c='orange',alpha=0.5)
ax.plot(size_check,dvz_check,label='dvz',c='b')
ax.plot(size_check,vturb_check,label='vturb',c='r')

ax.set_xlim(1e-5,1e3)
ax.set_ylim(1e-2,1e4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('a_p (cm)')
ax.set_ylabel('v_pp (cm/s)')
plt.legend()
plt.show()

srcterm=readcol(name+"/srcterm_inner.txt");
srcterm1=readcol(name+"/srcterm1au.txt");
srcterm3=readcol(name+"/srcterm3au.txt");
srcterm10=readcol(name+"/srcterm10au.txt");
srcterm30=readcol(name+"/srcterm30au.txt");
srcterm100=readcol(name+"/srcterm100au.txt");


m_peb_check=4./3.*np.pi*size_check**3*rho_peb




fig,ax=plt.subplots()

#ax.plot(size_check,srcterm)
ax.plot(m_peb_check,srcterm1*TUNIT,label='1 au')
ax.plot(m_peb_check,srcterm3*TUNIT,label='3 au')
ax.plot(m_peb_check,srcterm10*TUNIT,label='10 au')
ax.plot(m_peb_check,srcterm30*TUNIT,label='30 au')
ax.plot(m_peb_check,srcterm100*TUNIT,label='100 au')

ax.set_ylabel('srcterm*yr (g)')
ax.set_xlabel('grain mass (g)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.axis('equal')
plt.legend()
plt.show()


t001=np.linspace(1,1e4,1000001)
t01=np.linspace(1,1e4,100001)
t1=np.linspace(1,2e5,200001)
t10=np.linspace(1,2e5,20001)
t100=np.linspace(1,2e5,2001)
t1000=np.linspace(1,2e5,201)




fwd_eu_001=readcol("fwd_eu_001yr.txt")
fwd_eu_01=readcol("fwd_eu_01yr.txt")
#fwd_eu_1=readcol("fwd_eu_1yr.txt")
fwd_eu_10=readcol("fwd_eu_10yr.txt")
fwd_eu_100=readcol("fwd_eu_100yr.txt")
fwd_eu_1000=readcol("fwd_eu_1000yr.txt")




#rk2_01=readcol("rk2_01yr.txt")
#rk2_1=readcol("rk2_1yr.txt")
rk2_10=readcol("rk2_10yr.txt")
rk2_100=readcol("rk2_100yr.txt")
rk2_1000=readcol("rk2_1000yr.txt")


fig,ax=plt.subplots()




ax.plot(t001,fwd_eu_001,label="fwd euler, 0.01yr",c='C5')
ax.plot(t01,fwd_eu_01,label="fwd euler, 0.1yr",c='C4')
#ax.plot(t1,fwd_eu_1,label="fwd euler, 1yr",c='C0')
ax.plot(t10,fwd_eu_10,label="fwd euler, 10yr",c='C1')
ax.plot(t100,fwd_eu_100,label="fwd euler, 100yr",c='C2')
ax.plot(t1000,fwd_eu_1000,label="fwd euler, 1000yr",c='C3')




#ax.plot(t01,rk2_01,label="rk2, 0.1yr",ls='--',c='C4')
#ax.plot(t1,rk2_1,label="rk2, 1yr",ls='--',c='C0')
#ax.plot(t10,rk2_10,label="rk2, 10yr",ls='--',c='C1')
#ax.plot(t100,rk2_100,label="rk2, 100yr",ls='--',c='C2')
#ax.plot(t1000,rk2_1000,label="rk2, 1000yr",ls='--',c='C3')



ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()

plt.show()


tgrowth100au=readcol(name+"/tgrowth100au.txt");
tgrowth1au=readcol(name+"/tgrowth_inner.txt");

fig,ax=plt.subplots()

ax.plot(size_check,tgrowth100au/TUNIT,label='100au',c='b')
ax.plot(size_check,tgrowth1au/TUNIT,label='inner',c='r')

ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()

"""
tgrowth1au=readcol(name+"/tgrowth1au.txt");
tgrowth2au=readcol(name+"/tgrowth2au.txt");
tgrowth3au=readcol(name+"/tgrowth3au.txt");
tgrowth4au=readcol(name+"/tgrowth4au.txt");
tgrowth5au=readcol(name+"/tgrowth5au.txt");






fig,ax=plt.subplots()
ax.plot(size_check,tgrowth1au/TUNIT,label='1au')
ax.plot(size_check,tgrowth2au/TUNIT,label='2au')
ax.plot(size_check,tgrowth3au/TUNIT,label='3au')
ax.plot(size_check,tgrowth4au/TUNIT,label='4au')
ax.plot(size_check,tgrowth5au/TUNIT,label='5au')


ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()

vr1au=readcol(name+"/vr1au.txt");

fig,ax=plt.subplots()
ax.plot(size_check,-vr1au,label='1au')


#ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()
"""
gassig=readcol(name+"/gassig.txt")
gashei=readcol(name+"/gashei.txt")
gasrho=readcol(name+"/gasrho.txt")
yeta=readcol(name+"/yeta.txt")
yetavk=readcol(name+"/yetavk.txt")
cs=readcol(name+"/cs.txt")



fig=plt.figure(figsize=(6,7))
gs=gridspec.GridSpec(6,1)

ax1=plt.subplot(gs[0,0])
ax1.plot(r,gassig)
ax1.plot(r,1700*r**-1.5,ls='--')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel("gassig")


ax2=plt.subplot(gs[1,0])
ax2.plot(r,gashei/LUNIT)
ax2.plot(r,0.026*r**1.25,ls='--')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel("h_g")

ax3=plt.subplot(gs[2,0])
ax3.plot(r,gasrho)
ax3.plot(r,1.7e-9*r**(-11./4.),ls='--')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_ylabel("gasrho")

ax4=plt.subplot(gs[3,0])
ax4.plot(r,yeta)
ax4.plot(r,2.2e-3*r**0.5,ls='--')
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_ylabel("yeta")

ax5=plt.subplot(gs[4,0])
ax5.plot(r,yetavk)
ax5.plot(r,6600*r**0,ls='--')
#ax5.set_yscale('log')
ax5.set_xscale('log')
ax5.set_ylabel("yetavk")

ax6=plt.subplot(gs[5,0])
ax6.plot(r,cs)
ax6.plot(r,7.8e4*r**-0.25,ls='--')
ax6.set_yscale('log')
ax6.set_xscale('log')
ax6.set_ylabel("cs")

plt.show()
