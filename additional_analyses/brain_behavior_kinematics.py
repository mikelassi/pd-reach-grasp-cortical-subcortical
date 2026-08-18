"""
brain_behavior_kinematics.py
=========================================================================
Within-subject brain-behavior coupling using the FULL kinematic metric set.

Merges per-trial STN neural metrics (trial_table.csv) with the validated
reach-to-grasp kinematics (kinematic_metrics.csv) on subject/block/trial, then:
  per-subject Spearman rho(neural, behavior)  ->  group second-level sign-flip
  on Fisher-z (exhaustive 2^8). Behavior families pre-specified:
    PRIMARY (vigor/bradykinesia): peak velocity, movement duration.
    SECONDARY (smoothness/decomposition): submovement count, normalized jerk,
    smoothness index, amplitude.
=========================================================================
"""
import os, sys, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass
TEMP=os.path.dirname(os.path.abspath(__file__))
OUT=r"C:\Users\m.lassi\Documents\GitHub\pd-reach-grasp-cortical-subcortical\spectral_statistical_analysis\CBPT_results"

neu=pd.read_csv(os.path.join(TEMP,"trial_table.csv")).rename(columns={'subj':'subject'})
kin=pd.read_csv(r"F:\Projects\Parkinson_ReachGrasp\Reprocessing\kinematic_metrics.csv")
df=pd.merge(neu, kin, on=['subject','block','trial'], how='inner', suffixes=('','_k'))
df.to_csv(os.path.join(OUT,"brain_behavior_merged.csv"),index=False)
subs=sorted(df.subject.unique()); S=len(subs)
print(f"merged {len(df)} trials, {S} subjects -> brain_behavior_merged.csv")

def sf(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]; n=len(x)
    if n<4: return np.nan,np.nan
    sgn=np.array(np.meshgrid(*([[1,-1]]*n))).T.reshape(-1,n).astype(float)
    o=x.mean()/(x.std(ddof=1)/np.sqrt(n)+1e-12)
    pr=(sgn*x[None]).mean(1)/((sgn*x[None]).std(1,ddof=1)/np.sqrt(n)+1e-12)
    return o,np.mean(np.abs(pr)>=abs(o)-1e-9)
def subj_rhos(ncol,bcol):
    rs=[]
    for s in subs:
        d=df[df.subject==s][[ncol,bcol]].dropna()
        if len(d)>=6 and d[ncol].nunique()>2 and d[bcol].nunique()>2:
            r,_=st.spearmanr(d[ncol],d[bcol])
            if np.isfinite(r): rs.append(r)
    return np.array(rs)
def bh(p):
    p=np.asarray(p); n=p.size; o=np.argsort(p); q=np.empty(n); prev=1
    for r in range(n-1,-1,-1): prev=min(prev,p[o[r]]*n/(r+1)); q[o[r]]=prev
    return np.minimum(q,1)

NEURAL={'reach beta-ERD':'reach_beta_db','reach high-beta ERD':'reach_highbeta_db'}
# behavior: (label, column, family, expected-sign note)
BEHAV=[('peak velocity','reach_peak_vel_m_s','PRIMARY'),
       ('reach duration','dur_reach_s','PRIMARY'),
       ('submovements (n)','n_mvmt_reach','secondary'),
       ('normalized jerk','NJ_reach','secondary'),
       ('smoothness index','SI_reach','secondary'),
       ('reach amplitude','reach_median_radius_m','secondary')]

print("\nWITHIN-SUBJECT brain-behavior (per-subject Spearman -> group sign-flip on Fisher-z):")
rows=[]
for nlab,ncol in NEURAL.items():
    print(f"\n  [{nlab}]")
    ps=[]; recs=[]
    for blab,bcol,fam in BEHAV:
        rs=subj_rhos(ncol,bcol)
        z=np.arctanh(np.clip(rs,-.999,.999)); t,p=sf(z)
        recs.append((blab,fam,np.nanmean(rs) if len(rs) else np.nan,p,(rs<0).sum(),len(rs))); ps.append(p)
    q=bh(np.nan_to_num(np.array(ps),nan=1.0))
    for (blab,fam,mr,p,neg,n),qq in zip(recs,q):
        star='**' if p<0.05 else ''
        print(f"    {blab:18s} [{fam:9s}]: mean rho={mr:+.3f}  p={p:.4f} (FDR q={qq:.3f}){star}  {neg}/{n} neg")
        rows.append((nlab,blab,fam,mr,p,qq))

# figure: within-subject rho distributions for reach beta-ERD
fig,ax=plt.subplots(figsize=(11,5))
ncol='reach_beta_db'; data=[]; labs=[]
for blab,bcol,fam in BEHAV:
    rs=subj_rhos(ncol,bcol)
    data.append(list(rs)); labs.append(f"{blab}\n[{fam}]\nn={len(rs)}")
bp=ax.boxplot(data,patch_artist=True,showfliers=False)
for i,b in enumerate(bp['boxes']): b.set(facecolor='#aed6f1' if 'PRIMARY' in labs[i] else '#f5b7b1')
for i,d in enumerate(data): ax.scatter(np.full(len(d),i+1)+np.random.uniform(-.09,.09,len(d)),d,s=22,color='#1b4f72',zorder=3)
ax.axhline(0,color='k',lw=.8)
ax.set_xticklabels(labs,fontsize=8); ax.set_ylabel("within-subject Spearman rho")
ax.set_title("Reach STN beta-ERD vs trial kinematics (per subject; blue=primary, red=secondary)",fontweight='bold',fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"brain_behavior_kinematics.png"),dpi=180,bbox_inches='tight')
print(f"\nSaved {OUT}/brain_behavior_kinematics.png and brain_behavior_merged.csv")
