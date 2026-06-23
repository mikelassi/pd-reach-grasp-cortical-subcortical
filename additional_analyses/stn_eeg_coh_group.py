"""
stn_eeg_coh_group.py
=========================================================================
Group analysis of STN <-> scalp-EEG coherence by movement phase, with scalp
topographies and contralateral-vs-ipsilateral quantification.

CORRECTED inference: per-subject coherence (trials pooled within subject),
atanh-stabilised, group tested with exhaustive subject-level sign-flip on the
phase CONTRASTS (no trial-DOF significance).

LATERALITY (from A02_events_segmentation.m): wue05, wue09 moved the LEFT hand
(contralateral cortex = RIGHT hemisphere); all others moved the RIGHT hand.
Left-hand movers are mirror-flipped so that, for all subjects, the cortex
CONTRALATERAL to the moving hand is on the LEFT.

Input : stn_eeg_coh.mat (export_stn_eeg_coh.m): COH/ICOH [S x chan x freq x 3]
        conditions = {rest, reach, grip}
Output: CBPT_results/STN_EEG_coh_topo.png + console stats
=========================================================================
"""
import os, sys, re, warnings
import numpy as np, scipy.io as sio, scipy.stats as st
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import mne
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass

HERE=os.path.dirname(os.path.abspath(__file__))
OUT=os.path.join(os.path.dirname(HERE),"spectral_statistical_analysis","CBPT_results")
m=sio.loadmat(os.path.join(os.environ.get("TEMP",r"C:\Users\m.lassi\AppData\Local\Temp"),"stn_eeg_coh.mat"))
COH=np.asarray(m["COH"],float)          # [S, chan, freq, 3]
fr=m["freqs"].ravel().astype(float); labs=[str(x[0]).strip().lower() for x in m["labels"].ravel()]
subs=[str(x[0]) for x in m["SUBS"].ravel()]; S=COH.shape[0]
CONDS=["rest","reach","grip"]; LEFT={'wue05','wue09'}

# ---- mirror flip left-hand movers (contra cortex -> left) ----
def mirror(lab):
    mm=re.match(r'^([a-z]+)(\d+)(h?)$',lab)
    if not mm: return lab
    pre,num,h=mm.group(1),int(mm.group(2)),mm.group(3)
    if pre.endswith('z'): return lab
    return f"{pre}{num+1 if num%2==1 else num-1}{h}"
idx_of={l:i for i,l in enumerate(labs)}; mir=[idx_of.get(mirror(l),i) for i,l in enumerate(labs)]
COHf=COH.copy()
for s in range(S):
    if subs[s] in LEFT: COHf[s]=COH[s][mir,:,:]

# ---- beta band, atanh-stabilised ----
bm=(fr>=13)&(fr<=30)
beta=np.arctanh(np.sqrt(np.clip(COHf[:,:,bm,:].mean(2),0,0.999)))   # [S, chan, 3]  (atanh of |coh|)

info=mne.create_info(labs,400.,"eeg")
info.set_montage(mne.channels.make_standard_montage("standard_1005"),on_missing="ignore",match_case=False)
pos=np.array([info["chs"][i]["loc"][:3] for i in range(len(labs))])
valid=np.where(np.isfinite(pos).all(1)&(np.abs(pos).sum(1)>0))[0]
info2=mne.pick_info(info,valid)

sg=np.array(np.meshgrid(*([[1,-1]]*S))).T.reshape(-1,S).astype(float)
def sf(x,tail=0):
    o=x.mean()/(x.std(ddof=1)/np.sqrt(S)+1e-12); pr=(sg*x[None]).mean(1)/((sg*x[None]).std(1,ddof=1)/np.sqrt(S)+1e-12)
    return o,(np.mean(np.abs(pr)>=abs(o)-1e-9) if tail==0 else (np.mean(pr<=o+1e-9) if tail<0 else np.mean(pr>=o-1e-9)))

# ---- contra (left) vs ipsi (right) sensorimotor ROI ----
ROI_C=['c3','c1','c5','fcc3h','fcc1h','fcc5h','ccp3h','ccp1h','ccp5h','fc3','fc1','cp3','cp1']
ROI_I=[mirror(l) for l in ROI_C]
cidx=[i for i,l in enumerate(labs) if l in ROI_C]; iidx=[i for i,l in enumerate(labs) if l in ROI_I]
print("="*64); print("STN<->EEG beta coherence (n=8, contra=left after flip)"); print("="*64)
print(f"\nContralateral vs ipsilateral SM ROI (atanh|coh|), by phase:")
for cc,cn in enumerate(CONDS):
    c=beta[:,cidx,cc].mean(1); i=beta[:,iidx,cc].mean(1)
    o,p=sf(c-i,0)
    print(f"  {cn:6s}: contra={np.tanh(c.mean())**2:.3f} ipsi={np.tanh(i.mean())**2:.3f}  contra-ipsi t={o:+.2f} p={p:.4f}")
print("\nPhase modulation of CONTRALATERAL STN-cortex beta coherence (sign-flip):")
for a,b in [("reach","rest"),("grip","rest"),("reach","grip")]:
    x=beta[:,cidx,CONDS.index(a)].mean(1)-beta[:,cidx,CONDS.index(b)].mean(1)
    o,p=sf(x,0); print(f"  {a} - {b}: t={o:+.2f} p={p:.4f} ({int((x<0).sum())}/{S} down)")

# ---- topographies ----
def gm(cc): return np.tanh(beta[:,valid,cc].mean(0))**2          # back to coherence units
def gmdiff(a,b): return np.tanh(beta[:,valid,a].mean(0))**2 - np.tanh(beta[:,valid,b].mean(0))**2
fig,axes=plt.subplots(1,5,figsize=(19,4.2))
maps=[("rest",gm(0),"Reds",None),("reach",gm(1),"Reds",None),("grip",gm(2),"Reds",None),
      ("reach - rest",gmdiff(1,0),"RdBu_r",None),("grip - rest",gmdiff(2,0),"RdBu_r",None)]
for ax,(ttl,data,cmap,_) in zip(axes,maps):
    if cmap=="Reds": vlim=(0,np.nanmax(data))
    else: v=np.nanmax(np.abs(data)); vlim=(-v,v)
    im,_=mne.viz.plot_topomap(data,info2,axes=ax,cmap=cmap,vlim=vlim,contours=3,show=False)
    ax.set_title(f"STN-cortex beta coh\n{ttl}",fontsize=10)
    plt.colorbar(im,ax=ax,fraction=0.046,shrink=0.7)
fig.suptitle("STN <-> scalp-EEG beta coherence by reach-to-grasp phase (n=8; contralateral cortex = LEFT)",
             fontsize=13,fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT,"STN_EEG_coh_topo.png"),dpi=200,bbox_inches="tight")
print(f"\nSaved {OUT}/STN_EEG_coh_topo.png")
