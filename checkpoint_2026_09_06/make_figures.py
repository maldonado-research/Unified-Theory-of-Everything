"""Scientific illustrations from saved TOE diagnostic data; no fitted observations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
    'axes.spines.right':False,'axes.labelcolor':'#182b3a','text.color':'#182b3a',
    'axes.titleweight':'bold','figure.facecolor':'white','savefig.facecolor':'white'})
fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
s=np.logspace(-5,0,220)
G=2/3-2*s+2*s**1.5*np.arctan(1/np.sqrt(s))
pole=(2/9)/(s+1/3)
ax=axs[0]
ax.semilogx(s,G,lw=2.5,color='#006f85',label='Continuous spectrum')
ax.semilogx(s,pole,lw=2.2,ls='--',color='#bd6225',label='Matched one-pole model')
ax.set(xlabel='Euclidean squared-momentum variable s (toy units)',ylabel='Response G(s)',
       title='Same two moments; different response')
ax.legend(loc='lower left',frameon=False,fontsize=10)
ax.grid(alpha=.18)
ax.text(.04,.27,'Both match G(0) = 2/3 and G\N{PRIME}(0) = −2.\nOnly the continuum contains the s³ᐟ² term.',
        transform=ax.transAxes,fontsize=10,va='top')

data=json.loads((ROOT/'fermion_portal/BENCHMARK_RESULTS.json').read_text())
weak=next(c for c in data['cases'] if c['name']=='weak_certified')
k=np.array(weak['k_over_M']);n=np.array(weak['occupation_per_helicity'])
sel=(k>0)&(k<=3)
ax=axs[1]
ax.semilogy(k[sel],n[sel],lw=2.5,color='#006f85',label='Computed occupation per helicity')
ax.scatter([1],[1.25e-6],marker='^',color='#bd6225',s=55,zorder=5,label='Proved lower bound at k/M = 1')
ax.set(xlabel='Momentum k / asymptotic mass M',ylabel='Final occupation n(k)',
       title='Zero neutrino kernel; nonzero production',ylim=(1e-11,2e-5),xlim=(0,3))
ax.legend(loc='lower left',frameon=False,fontsize=9)
ax.grid(alpha=.18)
ax.text(.46,.93,'Prescribed weak mass pulse\na = 0.01, Mτ = 1\nNo mass crossing; no lepton asymmetry',
        transform=ax.transAxes,fontsize=9,va='top')
fig.suptitle('TOE applicability checks: spectral information and time dependence',fontsize=16,fontweight='bold')
fig.supxlabel('Illustrative models and analytic controls. No experimental data or calibrated cosmological prediction.',fontsize=10)
fig.savefig(ROOT/'TOE_APPLICABILITY_CHECKS.png',dpi=180)
fig.savefig(ROOT/'TOE_APPLICABILITY_CHECKS.svg')
plt.close(fig)

fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
colors=['#006f85','#bd6225','#644d92','#29815b']
names=['shallow_fast','shallow_slow','crossing_fast','crossing_slow']
for name,color in zip(names,colors):
    c=next(c for c in data['cases'] if c['name']==name)
    kk=np.array(c['k_over_M']);nn=np.array(c['occupation_per_helicity'])
    sel=(kk>0)&(kk<=4)
    ax.semilogy(kk[sel],nn[sel],lw=2,color=color,label=f"a = {c['pulse_depth_a']:g}, Mτ = {c['width_M_tau']:g}")
ax.set(title='Pulse history changes the heavy-fermion spectrum',xlabel='Momentum k / M',
       ylabel='Final occupation per helicity',ylim=(1e-10,1.2),xlim=(0,4))
ax.legend(frameon=False,ncol=2)
ax.grid(alpha=.18)
fig.supxlabel('Flat-spacetime toy; same positive asymptotic mass. No backreaction or cosmological fit.',fontsize=10)
fig.savefig(ROOT/'fermion_portal/PULSE_SPECTRA.png',dpi=180)
fig.savefig(ROOT/'fermion_portal/PULSE_SPECTRA.svg')
