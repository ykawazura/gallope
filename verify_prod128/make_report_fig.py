import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

INJ=1.0; nexp_m=6; nexp_p=4; nu_m=1.0; mu_p=10.0

# --- full W_free(t) trajectory (all 4 segments, t=0.5 -> 21.76) ---
segs=['echo128_seg_2355202_t4.nc','echo128_seg_2357102_t4to10.nc',
      'echo128_seg_2358620.nc','echo128_seg_2358914.nc']
T=[];W=[]
for f in segs:
    d=Dataset(f); v=d.variables
    tt=np.asarray(v['tt'][:]); wf=np.asarray(v['W_free'][:])
    m=(tt>1e-6)|(wf>1e-6)   # drop spurious t=0 zero rows
    T.append(tt[m]);W.append(wf[m]); d.close()
T=np.concatenate(T);W=np.concatenate(W); o=np.argsort(T);T=T[o];W=W[o]
# keep physical branch (drop the raw-IC t=0 spike)
keep=T>=0.5; T=T[keep];W=W[keep]

# --- flux + spectra from the STEADY-window segment (t=16.5->21.76, flux ON) ---
d=Dataset('echo128_seg_2358914.nc'); v=d.variables
mm=np.asarray(v['mm'][:]); kp=np.asarray(v['kpbin'][:])
Gk=np.asarray(v['Gamma_m_kint'][:]); g2=np.asarray(v['g2_bin'][:]) # (tt,kz,mm,kpbin)
Wm=np.asarray(v['W_m'][:]) if 'W_m' in v else None
n=Gk.shape[0]; sl=slice(n//2,n)
Gt=Gk[sl].mean(axis=0)
Wmt=Wm[sl].mean(axis=0) if Wm is not None else g2[sl].mean(axis=(0,1,3))
# Elsasser perp spectrum (turbulence health)
zpp=np.asarray(v['zppe2_bin'][:])[sl].mean(axis=(0,1)) if 'zppe2_bin' in v else None
zmp=np.asarray(v['zmpe2_bin'][:])[sl].mean(axis=(0,1)) if 'zmpe2_bin' in v else None
d.close()

nm=mm.size
fig,ax=plt.subplots(2,2,figsize=(12,8.2))

# (A) k-int Hermite flux = THE result
a=ax[0,0]
a.axhline(INJ,ls='--',c='crimson',lw=1,label='injection rate = 1')
a.axhline(0,ls=':',c='k',lw=1,label='echo / fluidized (=0)')
a.axvspan(48,nm,color='0.85',label=r'$\nu m^6$ dissipation range')
a.plot(mm,Gt,'-',c='C0',lw=2,label=r'$\langle\Gamma(m)\rangle_t$  (k$_\perp$+k$_z$ integrated)')
a.set_xlim(0,nm); a.set_ylim(-0.15,1.2)
a.set_xlabel('m (Hermite)'); a.set_ylabel(r'$\Gamma(m)$')
a.set_title('(A) Hermite flux: FORWARD cascade (phase mixing), NOT echo',fontsize=10,weight='bold')
a.legend(fontsize=7,loc='lower left')
a.annotate('plateau $\\approx$ injection\n(=0.98 across m=20-48)',(30,1.05),fontsize=8,color='C0')
a.annotate('telescopes to 0\nat top m',(95,0.12),fontsize=7,color='0.3')

# (B) Hermite spectrum W_m
b=ax[0,1]
b.loglog(mm[1:],Wmt[1:],'o-',c='C0',ms=3,lw=1,label=r'$\langle W_m\rangle$')
mr=mm[(mm>=2)&(mm<=48)]
b.loglog(mr,Wmt[2]* (mr/2.0)**-0.5,'--',c='crimson',lw=1,label=r'$m^{-1/2}$ (phase-mix)')
b.loglog(mr,Wmt[2]* (mr/2.0)**-1.0,'--',c='green',lw=1,label=r'$m^{-1}$ (fluidized)')
b.set_xlabel('m'); b.set_ylabel(r'$W_m$')
b.set_title('(B) Hermite free-energy spectrum',fontsize=10,weight='bold')
b.legend(fontsize=8)

# (C) W_free(t) full driven trajectory -> steady state
c=ax[1,0]
c.plot(T,W,'o-',c='C0',ms=3,lw=1)
c.axvspan(0.5,2.5,color='0.9')
# steady window used for the A/B/D averages (later half of the t=16.5-21.76 segment)
steady=(T>=19.0)
c.axvspan(19.0,T[-1],color='#cfe8d4')
Wsteady=W[steady].mean()
c.axhline(Wsteady,ls='--',c='green',lw=1)
c.annotate('IC transient',(1.5,W.max()*0.999),fontsize=7,ha='center',color='0.4')
c.annotate('driven rise\n(escaped free decay)',(8,W.min()+0.3),fontsize=8,color='C0')
c.annotate('peak $\\approx$17.66 (t$\\approx$19.5)\nthen settles: plateau\n%.2f$\\pm$%.2f (0.3%%)'
           %(Wsteady,W[steady].std()),(15.0,W.min()+0.15),fontsize=7.5,color='green')
c.set_xlabel('t'); c.set_ylabel(r'$W_{free}$')
c.set_title('(C) Free energy: reaches statistical steady state (Gate 2)',fontsize=10,weight='bold')

# (D) Elsasser perp spectrum
dd=ax[1,1]
if zpp is not None:
    msk=kp>=1
    dd.loglog(kp[msk],zpp[msk],'o-',c='C3',ms=3,lw=1,label=r'$|z^+|^2$')
    dd.loglog(kp[msk],zmp[msk],'s-',c='C0',ms=3,lw=1,label=r'$|z^-|^2$')
    kref=kp[(kp>=3)&(kp<=20)]
    dd.loglog(kref, zpp[msk][ (kp[msk]>=3)&(kp[msk]<=20) ][0]*(kref/kref[0])**(-5/3.),
              '--',c='0.4',lw=1,label=r'$k_\perp^{-5/3}$')
    dd.legend(fontsize=8)
dd.set_xlabel(r'$k_\perp$'); dd.set_ylabel('perp spectrum')
dd.set_title('(D) Alfvenic (Elsasser) perp spectrum: turbulence health',fontsize=10,weight='bold')

fig.suptitle('prod128 driven KRMHD (128$^2$x64, nm=128, ene_inj_g=1): phase-mixing-dominated forward Hermite cascade  [t=%.1f-%.1f]'%(T[0],T[-1]),
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig('report_prod128_fig.png',dpi=110)
print("saved report_prod128_fig.png  (t=%.2f..%.2f)"%(T[0],T[-1]))
print("Gamma plateau m20-48 = %.3f ; W_m slope (m=3..40) = %.3f"%(
    Gt[(mm>=20)&(mm<=48)].mean(),
    np.polyfit(np.log(mm[(mm>=3)&(mm<=40)]),np.log(Wmt[(mm>=3)&(mm<=40)]),1)[0]))
