"""
brain_behavior_clinic.py
=========================================================================
Reusable brain - behavior - clinic framework.

  WITHIN-SUBJECT brain-behavior (the POWERED analysis):
    For each subject, correlate trial-level neural state (reach beta ERD) with
    trial-level behavior (peak velocity, movement duration). Aggregate the
    per-subject correlations to a GROUP second-level test (Fisher-z of the
    within-subject rho, exhaustive subject-level sign-flip): is the
    brain-behavior coupling consistent across patients?

  BETWEEN-SUBJECT brain-clinic (EXPLORATORY, n=8):
    Per-subject neural / behavioral summaries vs clinical scores (UPDRS off/off,
    disease duration, DAT asymmetry). Reported as hypothesis-generating only.

Inputs : trial_table.csv (export_trial_table.m)  -- enrich with more kinematic
         columns when the detailed extraction syncs; the framework picks them up.
         ClinicalData.xlsx
Output : CBPT_results/brain_behavior_clinic.png + console
=========================================================================
"""
import os, sys, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass

TEMP=os.path.dirname(os.path.abspath(__file__))
OUT=r"C:\Users\m.lassi\Documents\GitHub\pd-reach-grasp-cortical-subcortical\spectral_statistical_analysis\CBPT_results"
T=pd.read_csv(os.path.join(TEMP,"trial_table.csv"))
subs=sorted(T.subj.unique())
S=len(subs); sg=np.array(np.meshgrid(*([[1,-1]]*S))).T.reshape(-1,S).astype(float)
def sf(x):
    x=np.asarray(x); o=x.mean()/(x.std(ddof=1)/np.sqrt(S)+1e-12)
    pr=(sg*x[None]).mean(1)/((sg*x[None]).std(1,ddof=1)/np.sqrt(S)+1e-12)
    return o,np.mean(np.abs(pr)>=abs(o)-1e-9)

# ---------- WITHIN-SUBJECT brain-behavior ----------
NEURAL=["reach_beta_db","reach_highbeta_db"]
BEHAV =["peak_reach_vel","mvt_dur_s","reach_dur_s"]
print("="*70); print("WITHIN-SUBJECT brain-behavior (per-subject Spearman -> group sign-flip)")
print("hypothesis: deeper reach beta-ERD -> faster / shorter movement"); print("="*70)
results={}
for nf in NEURAL:
    for bf in BEHAV:
        rs=[]
        for s in subs:
            d=T[T.subj==s][[nf,bf]].dropna()
            if len(d)>=6:
                r,_=st.spearmanr(d[nf],d[bf]); rs.append(r)
        rs=np.array(rs); z=np.arctanh(np.clip(rs,-0.999,0.999))
        t,p=sf(z); results[(nf,bf)]=(rs,p)
        print(f"  {nf:18s} x {bf:15s}: mean rho={rs.mean():+.3f}  group p={p:.4f}  ({np.sum(rs<0)}/{len(rs)} neg)")

# ---------- BETWEEN-SUBJECT brain-clinic ----------
df=pd.read_excel(r"F:\Projects\Parkinson_ReachGrasp\Raw\ClinicalData.xlsx").rename(columns={'Unnamed: 0':'id'})
cl=df[df['id'].astype(str).str.startswith('wue')].set_index('id')
def cnum(s,col):
    try: return float(cl.loc[s,col])
    except: return np.nan
clin=pd.DataFrame({s:dict(
    UPDRS=cnum(s,'UPDRS-III post-DBS meds-off, stim-off (score)'),
    duration=cnum(s,'Disease duration at surgery (years)'),
    DATasym=cnum(s,'Unnamed: 16')) for s in subs}).T
# per-subject neural/behavioral summaries
summ=T.groupby('subj').agg(reach_beta=('reach_beta_db','mean'),
                           mvt_dur=('mvt_dur_s','median'),
                           peak_vel=('peak_reach_vel','median')).loc[subs]
M=summ.join(clin)
print("\n"+"="*70); print("BETWEEN-SUBJECT brain/behav-clinic (Spearman, n=8, EXPLORATORY)"); print("="*70)
clin_corr=[]
for nf in ["reach_beta","mvt_dur","peak_vel"]:
    for cf in ["UPDRS","duration","DATasym"]:
        ok=M[[nf,cf]].dropna()
        r,p=st.spearmanr(ok[nf],ok[cf]) if len(ok)>=4 else (np.nan,np.nan)
        clin_corr.append((nf,cf,r,p,len(ok)))
        print(f"  {nf:12s} x {cf:10s}: rho={r:+.2f} p={p:.3f} (n={len(ok)})")

# ---------- figure ----------
fig=plt.figure(figsize=(14,7));
# panel A: within-subject rho distributions
ax=fig.add_subplot(2,3,1)
labels=[]; data=[]
for (nf,bf),(rs,p) in results.items():
    labels.append(f"{nf.split('_')[1]}\nx {bf.split('_')[0]}"); data.append(rs)
bp=ax.boxplot(data,patch_artist=True,showfliers=False)
for b in bp['boxes']: b.set(facecolor="#aed6f1")
for i,d in enumerate(data): ax.scatter(np.full(len(d),i+1)+np.random.uniform(-.1,.1,len(d)),d,s=18,color="#21618c",zorder=3)
ax.axhline(0,color='k',lw=.8); ax.set_xticklabels(labels,fontsize=7,rotation=20)
ax.set_ylabel("within-subject Spearman rho"); ax.set_title("Brain-behavior (per subject)",fontsize=10,fontweight='bold')
# panel B: example within-subject scatter (best pair)
best=min(results.items(),key=lambda kv: kv[1][1]); (nf,bf)=best[0]
ax=fig.add_subplot(2,3,2)
for s in subs:
    d=T[T.subj==s][[nf,bf]].dropna()
    if len(d)>=6: ax.scatter(d[nf],d[bf],s=10,alpha=.5,label=s[3:])
ax.set_xlabel(nf); ax.set_ylabel(bf); ax.set_title(f"trials: {nf} x {bf}\ngroup p={best[1][1]:.3f}",fontsize=9)
# panel C-F: clinic scatters
for k,(nf,cf,r,p,n) in enumerate([c for c in clin_corr][:4]):
    ax=fig.add_subplot(2,3,3+k)
    ax.scatter(M[nf],M[cf],s=55,color="#8e44ad",edgecolor='k')
    for s in subs:
        if np.isfinite(M.loc[s,nf]) and np.isfinite(M.loc[s,cf]): ax.annotate(s[3:],(M.loc[s,nf],M.loc[s,cf]),fontsize=7)
    ax.set_xlabel(nf,fontsize=8); ax.set_ylabel(cf,fontsize=8)
    ax.set_title(f"rho={r:+.2f} p={p:.3f} (n={n})",fontsize=9,color=('#1e8449' if p<0.05 else '#666'))
fig.suptitle("Brain-Behavior-Clinic framework  (within-subject behavior = powered; clinic = exploratory n=8)",
             fontsize=12,fontweight='bold')
fig.tight_layout()
fig.savefig(os.path.join(OUT,"brain_behavior_clinic.png"),dpi=170,bbox_inches='tight')
print(f"\nSaved {OUT}/brain_behavior_clinic.png")
