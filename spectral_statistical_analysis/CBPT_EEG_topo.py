"""
CBPT_EEG_topo.py
=========================================================================
EEG topographic analysis of the reach-to-grasp task (analogous to the LFP
spectral statistics), using a CLEAN pre-movement baseline.

  (a) baseline = early part of the pre-first-movement rest (movement-prep
      window trimmed), to avoid normalising away the anticipatory ERD.
  (b) per-channel band power -> group spatial cluster test (sign-flip, channel
      adjacency from the montage) -> scalp topographies of where the reach ERD
      and the post-movement beta rebound (PMBR) sit.

Input : eeg_topo.mat  (export_eeg_topo.m) — reach/rest/restpost PSD [S x chan x freq]
Output: CBPT_results/EEG_topo.png + console cluster stats

Author: Michael Lassi
=========================================================================
"""
import os, sys, warnings
import numpy as np, scipy.io as sio, scipy.stats as st
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import mne
from mne.stats import permutation_cluster_1samp_test
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"CBPT_results")
m=sio.loadmat(os.path.join(HERE,"eeg_topo.mat"))
fr=m["freqs"].ravel().astype(float); labs=[str(x[0]).strip() for x in m["labels"].ravel()]
R=np.asarray(m["reach_psd"],float); B=np.asarray(m["rest_psd"],float); Pp=np.asarray(m["restpost_psd"],float)
S=R.shape[0]

info=mne.create_info(labs,400.,"eeg")
info.set_montage(mne.channels.make_standard_montage("standard_1005"),on_missing="ignore",match_case=False)
pos=np.array([info["chs"][i]["loc"][:3] for i in range(len(labs))])
valid=np.where(np.isfinite(pos).all(1)&(np.abs(pos).sum(1)>0))[0]
info2=mne.pick_info(info,valid)
adj,_=mne.channels.find_ch_adjacency(info2,"eeg")
print(f"{len(valid)}/{len(labs)} channels with positions")

def cluster(NUM,DEN,lo,hi,tail=0):
    fm=(fr>=lo)&(fr<=hi)
    X=np.array([10*np.log10(NUM[s][np.ix_(valid,fm)].mean(1)/DEN[s][np.ix_(valid,fm)].mean(1)) for s in range(S)])
    thr=st.t.ppf(0.975,S-1)
    T,cl,pv,_=permutation_cluster_1samp_test(X,threshold=thr,n_permutations=2**S,tail=tail,
                                             adjacency=adj,out_type="mask",seed=0,verbose=False)
    mask=np.zeros(len(valid),bool); sig=[]
    for c,p in zip(cl,pv):
        if p<0.05: mask|=c; sig.append((int(c.sum()),float(p)))
    minp=float(min(pv)) if len(pv) else 1.0
    return X.mean(0),mask,sig,minp

PANELS=[("Reach vs rest","mu (8-13 Hz)",R,B,8,13),
        ("Reach vs rest","beta (13-30 Hz)",R,B,13,30),
        ("PMBR (post vs rest)","beta (13-30 Hz)",Pp,B,13,30),
        ("PMBR (post vs rest)","high-beta (20-30 Hz)",Pp,B,20,30)]

fig,axes=plt.subplots(1,4,figsize=(16,4.4))
print("\nTopographic spatial cluster tests (n=8):")
for ax,(name,bn,NUM,DEN,lo,hi) in zip(axes,PANELS):
    gm,mask,sig,minp=cluster(NUM,DEN,lo,hi)
    vlim=np.nanmax(np.abs(gm))
    im,_=mne.viz.plot_topomap(gm,info2,axes=ax,cmap="RdBu_r",vlim=(-vlim,vlim),
                              mask=mask,mask_params=dict(marker='o',markerfacecolor='k',
                              markeredgecolor='k',markersize=4),contours=4,show=False)
    sigtxt=("; ".join(f"p={p:.3f} ({n}ch)" for n,p in sig)) if sig else f"n.s. (min p={minp:.2f})"
    ax.set_title(f"{name}\n{bn}\n{sigtxt}",fontsize=10)
    plt.colorbar(im,ax=ax,fraction=0.046,shrink=0.7,label="dB")
    print(f"  {name:20s} {bn:20s}: {sigtxt}")
fig.suptitle("EEG scalp topography (n=8, clean pre-movement baseline) - black dots = significant cluster channels",
             fontsize=13,fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT,"EEG_topo.png"),dpi=200,bbox_inches="tight")
print(f"\nSaved {OUT}/EEG_topo.png")
