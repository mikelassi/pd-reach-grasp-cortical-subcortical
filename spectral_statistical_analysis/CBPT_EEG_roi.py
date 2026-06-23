"""
CBPT_EEG_roi.py
=========================================================================
ROI-restricted EEG topographic cluster test (more powerful than whole-scalp
at n=8), with hemisphere-flipping so the cortex CONTRALATERAL to the moving
hand is aligned across subjects.

  - wue05, wue09 moved the LEFT hand -> their channels are mirror-flipped L<->R
    so that, for all subjects, the contralateral sensorimotor cortex sits on
    the LEFT. (Right-hand movers are left as-is.)
  - Spatial sign-flip cluster test restricted to a left/midline sensorimotor
    ROI (~22 channels) instead of all 124.

Output: CBPT_results/EEG_roi.png + console stats
=========================================================================
"""
import os, sys, re, warnings
import numpy as np, scipy.io as sio, scipy.stats as st
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import mne
from mne.stats import permutation_cluster_1samp_test
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"CBPT_results")
m=sio.loadmat(os.path.join(HERE,"eeg_topo.mat"))
fr=m["freqs"].ravel().astype(float); labs=[str(x[0]).strip().lower() for x in m["labels"].ravel()]
R=np.asarray(m["reach_psd"],float); B=np.asarray(m["rest_psd"],float); Pp=np.asarray(m["restpost_psd"],float)
subs=[str(x[0]) for x in m["SUBJECTS"].ravel()]; S=R.shape[0]
LEFT={'wue05','wue09'}

# ---- left<->right mirror map ----
def mirror(lab):
    mm=re.match(r'^([a-z]+)(\d+)(h?)$',lab)
    if not mm: return lab
    pre,num,h=mm.group(1),int(mm.group(2)),mm.group(3)
    if pre.endswith('z'): return lab
    new=num+1 if num%2==1 else num-1
    return f"{pre}{new}{h}"
idx_of={l:i for i,l in enumerate(labs)}
mir_idx=[idx_of.get(mirror(l),i) for i,l in enumerate(labs)]   # index to copy from when flipping

def flip(arr):  # flip the left-hand movers
    out=arr.copy()
    for s in range(S):
        if subs[s] in LEFT: out[s]=arr[s][mir_idx,:]
    return out
Rf,Bf,Ppf=flip(R),flip(B),flip(Pp)

# ---- montage + ROI (contralateral = LEFT after flip) ----
info=mne.create_info(labs,400.,"eeg")
info.set_montage(mne.channels.make_standard_montage("standard_1005"),on_missing="ignore",match_case=False)
ROI=['fc5','fc3','fc1','fcz','c5','c3','c1','cz','cp5','cp3','cp1','cpz',
     'fcc5h','fcc3h','fcc1h','ccp5h','ccp3h','ccp1h','ffc1h','ffc3h','cpp1h','cpp3h']
pos=np.array([info["chs"][i]["loc"][:3] for i in range(len(labs))])
roi=[i for i,l in enumerate(labs) if l in ROI and np.isfinite(pos[i]).all() and np.abs(pos[i]).sum()>0]
info_roi=mne.pick_info(info,roi); roi_labs=[labs[i] for i in roi]
adj,_=mne.channels.find_ch_adjacency(info_roi,"eeg")
print(f"ROI channels: {len(roi)}")

sg=np.array(np.meshgrid(*([[1,-1]]*S))).T.reshape(-1,S).astype(float)
def roi_cluster(NUM,DEN,lo,hi,tail=0):
    fm=(fr>=lo)&(fr<=hi)
    X=np.array([10*np.log10(NUM[s][np.ix_(roi,fm)].mean(1)/DEN[s][np.ix_(roi,fm)].mean(1)) for s in range(S)])
    thr=st.t.ppf(0.975,S-1)
    T,cl,pv,_=permutation_cluster_1samp_test(X,threshold=thr,n_permutations=2**S,tail=tail,
                                             adjacency=adj,out_type="mask",seed=0,verbose=False)
    mask=np.zeros(len(roi),bool); sig=[]
    for c,p in zip(cl,pv):
        if p<0.05: mask|=c; sig.append((int(c.sum()),float(p)))
    # also ROI-mean test
    xm=X.mean(1); o=xm.mean()/(xm.std(ddof=1)/np.sqrt(S)+1e-12)
    pr=(sg*xm[None]).mean(1)/((sg*xm[None]).std(1,ddof=1)/np.sqrt(S)+1e-12)
    pmean=np.mean(np.abs(pr)>=abs(o)-1e-9)
    return X.mean(0),mask,sig,(float(min(pv)) if len(pv) else 1.0),xm.mean(),pmean

PANELS=[("Reach vs rest","mu (8-13 Hz)",Rf,Bf,8,13),
        ("Reach vs rest","beta (13-30 Hz)",Rf,Bf,13,30),
        ("PMBR (post vs rest)","beta (13-30 Hz)",Ppf,Bf,13,30),
        ("PMBR (post vs rest)","high-beta (20-30 Hz)",Ppf,Bf,20,30)]
fig,axes=plt.subplots(1,4,figsize=(16,4.6))
print("\nROI-restricted (hemisphere-flipped, contralateral=LEFT) cluster tests, n=8:")
for ax,(name,bn,NUM,DEN,lo,hi) in zip(axes,PANELS):
    gm,mask,sig,minp,roimean,pmean=roi_cluster(NUM,DEN,lo,hi)
    vlim=np.nanmax(np.abs(gm))
    im,_=mne.viz.plot_topomap(gm,info_roi,axes=ax,cmap="RdBu_r",vlim=(-vlim,vlim),
                              mask=mask,mask_params=dict(marker='o',markerfacecolor='k',markersize=6),
                              contours=0,sphere=0.11,show=False)
    cltxt=("; ".join(f"clust p={p:.3f}({n}ch)" for n,p in sig)) if sig else f"clust n.s.(min p={minp:.2f})"
    ax.set_title(f"{name}\n{bn}\nROI-mean {roimean:+.2f}dB p={pmean:.3f}\n{cltxt}",fontsize=9)
    print(f"  {name:20s} {bn:16s}: ROI-mean {roimean:+.2f}dB p={pmean:.3f} | {cltxt}")
fig.suptitle("EEG sensorimotor ROI (hemisphere-flipped, contralateral=left) — n=8",fontsize=13,fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT,"EEG_roi.png"),dpi=200,bbox_inches="tight")
print(f"\nSaved {OUT}/EEG_roi.png")
