"""
CBPT_LFP_subjectlevel.py
=========================================================================
Cluster-based permutation test (CBPT) on STN-LFP phase spectra — corrected,
subject-level version.

Fixes the three problems in CBPT_LFP_raw.ipynb:
  1. PSEUDOREPLICATION — the old test pooled ~239 trials across subjects and
     treated them as independent observations (df ~238). Here the statistical
     unit is the SUBJECT (n=8, df=7): one median spectrum per subject per
     phase, paired differences, one-sample sign-flip permutation. This is the
     Maris & Oostenveld (2007) design and the level at which a group claim is
     valid.
  2. DOUBLE-DIPPING — the old code selected bands with the CBPT, then re-tested
     those same bands with paired t-tests on the same data. Removed: the
     cluster-mass p-value IS the multiple-comparison-corrected inference. No
     follow-up tests.
  3. EDGE ARTIFACTS / axis drift — the old code interpolated+extrapolated each
     subject onto a uniform grid (fabricating power at 1 and 80 Hz). v2 already
     gives every subject an identical cwtfilterbank frequency axis, so we use it
     directly with no interpolation.

Inference details
-----------------
  * Cluster-forming threshold: two-tailed parametric t, df = n_subj-1.
  * Permutation: exhaustive sign-flip. With n=8 there are only 2^8 = 256
    sign permutations, so the SMALLEST achievable cluster p-value is
    1/256 ~= 0.0039. p-values below that are impossible under the correct null
    (the old test's p=0.0001 was a symptom of trial pooling).
  * Frequency adjacency: 1-D chain (neighbouring frequency bins), via
    adjacency=None (MNE treats the single remaining axis as a regular lattice).
  * Phases (from v2, 5 of them): Rest_pre is the normalisation baseline (~0 dB
    by construction); Rest_post carries the post-movement beta rebound (PMBR).
  * Family-wise correction across the reported contrasts is reported as a
    Bonferroni-adjusted cluster p (p_fwer) alongside the raw cluster p.

Inputs : cbpt_input.mat   (from export_psd_for_cbpt.m)
Outputs: CBPT_results/    PNG per contrast + overview + clusters.json + summary.txt

Author: Michael Lassi
=========================================================================
"""

import os
import json
import itertools
import numpy as np
import scipy.io as sio
import scipy.stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mne.stats import permutation_cluster_1samp_test

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
IN_MAT = os.path.join(HERE, "cbpt_input.mat")
OUT_DIR = os.path.join(HERE, "CBPT_results")
os.makedirs(OUT_DIR, exist_ok=True)

ALPHA_CLUSTER = 0.05   # cluster-forming threshold alpha (two-tailed)
ALPHA_SIG = 0.05       # cluster-level significance
FREQ_MAX_PLOT = 80     # Hz

PHASE_COLORS = {
    "Rest_pre":  "#c0392b",
    "Reach":     "#27ae60",
    "Grasp":     "#2980b9",
    "Pull":      "#8e44ad",
    "Rest_post": "#e67e22",
}

# A-priori contrast families (reported grouped). Each entry: (A, B, family)
CONTRASTS = [
    ("Rest_pre", "Reach",     "Movement vs baseline (ERD/ERS)"),
    ("Rest_pre", "Grasp",     "Movement vs baseline (ERD/ERS)"),
    ("Rest_pre", "Pull",      "Movement vs baseline (ERD/ERS)"),
    ("Rest_pre", "Rest_post", "Post-movement beta rebound (PMBR)"),
    ("Reach",    "Grasp",     "Within-movement progression"),
    ("Grasp",    "Pull",      "Within-movement progression"),
    ("Reach",    "Pull",      "Within-movement progression"),
]


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------
def _cellstr(a):
    out = []
    for x in np.asarray(a).ravel():
        if isinstance(x, np.ndarray):
            out.append(str(x.ravel()[0]) if x.size else "")
        else:
            out.append(str(x))
    return out


m = sio.loadmat(IN_MAT)
freqs = m["freqs"].ravel().astype(float)
subjects = _cellstr(m["subjects"])
phase_names = _cellstr(m["phase_names"])
psd = np.asarray(m["psd_subj_phase"], dtype=float)   # [n_subj, n_phase, n_fr]
n_trials_subj = m["n_trials_subj"].ravel().astype(int)
n_subj, n_phase, n_fr = psd.shape
pidx = {p: i for i, p in enumerate(phase_names)}

df = n_subj - 1
thresh = scipy.stats.t.ppf(1 - ALPHA_CLUSTER / 2, df)   # two-tailed
n_perm_exact = 2 ** n_subj
p_floor = 1.0 / n_perm_exact

print("=" * 68)
print("Subject-level CBPT on STN-LFP phase spectra")
print("=" * 68)
print(f"Subjects (n={n_subj}): {', '.join(subjects)}")
print(f"Trials pooled/subject : {n_trials_subj.tolist()}  (total {n_trials_subj.sum()})")
print(f"Frequencies           : {n_fr} bins, {freqs.min():.2f}-{freqs.max():.2f} Hz")
print(f"Cluster-forming t      : |t| > {thresh:.3f}  (df={df}, two-tailed a={ALPHA_CLUSTER})")
print(f"Permutation            : exhaustive sign-flip, {n_perm_exact} perms "
      f"(p floor = {p_floor:.4f})")
print("=" * 68)


# ----------------------------------------------------------------------
# Run CBPT per contrast
# ----------------------------------------------------------------------
def _to_mask(cl, n):
    """Normalise an MNE cluster into a 1-D boolean mask of length n.
    Handles boolean masks, integer index arrays, and tuples of slices/arrays
    (MNE's 1-D clusters come back as a 1-tuple holding a slice)."""
    mask = np.zeros(n, dtype=bool)
    el = cl[0] if isinstance(cl, tuple) else cl
    if isinstance(el, slice):
        mask[el] = True
    else:
        arr = np.asarray(el).ravel()
        if arr.dtype == bool:
            mask[: arr.size] = arr
        else:
            mask[arr.astype(int)] = True
    return mask


def run_cbpt(phaseA, phaseB):
    """One-sample sign-flip CBPT on subject-level paired differences A-B."""
    X = psd[:, pidx[phaseA], :] - psd[:, pidx[phaseB], :]   # [n_subj, n_fr]
    T_obs, clusters, cluster_pv, H0 = permutation_cluster_1samp_test(
        X,
        threshold=thresh,
        n_permutations=n_perm_exact,   # exact sign-flip for small n
        tail=0,
        adjacency=None,                # 1-D frequency chain
        out_type="mask",
        seed=0,
        verbose=False,
    )
    return X, T_obs, clusters, np.asarray(cluster_pv)


results = []
n_contrasts = len(CONTRASTS)

for (A, B, family) in CONTRASTS:
    X, T_obs, clusters, cluster_pv = run_cbpt(A, B)
    sig = []
    for cl, pv in zip(clusters, cluster_pv):
        mask = _to_mask(cl, n_fr)
        f_in = freqs[mask]
        rec = {
            "f_lo": float(f_in.min()),
            "f_hi": float(f_in.max()),
            "n_bins": int(mask.sum()),
            "mass": float(T_obs[mask].sum()),
            "t_peak": float(T_obs[mask][np.argmax(np.abs(T_obs[mask]))]),
            "p": float(pv),
            "p_fwer": float(min(pv * n_contrasts, 1.0)),
            "sign": "A>B" if T_obs[mask].sum() < 0 else "A<B",  # X=A-B; +T => A>B
            "mask": mask,
        }
        # X = A - B  => positive T means A > B
        rec["sign"] = "A>B" if rec["t_peak"] > 0 else "A<B"
        sig.append(rec)
    results.append({"A": A, "B": B, "family": family,
                    "X": X, "T_obs": T_obs, "clusters": sig})


# ----------------------------------------------------------------------
# Console report
# ----------------------------------------------------------------------
lines = []
def emit(s=""):
    print(s); lines.append(s)

emit("\nCLUSTER RESULTS  (sig = cluster p < %.2f)" % ALPHA_SIG)
emit("-" * 68)
cur_family = None
for r in results:
    if r["family"] != cur_family:
        cur_family = r["family"]
        emit(f"\n## {cur_family}")
    head = f"{r['A']} vs {r['B']}"
    sig_cl = [c for c in r["clusters"] if c["p"] < ALPHA_SIG]
    if not sig_cl:
        emit(f"  {head:24s}  no significant cluster")
        continue
    for c in sig_cl:
        direction = f"{r['A']}>{r['B']}" if c['t_peak'] > 0 else f"{r['A']}<{r['B']}"
        floor = "  [at p-floor]" if abs(c["p"] - p_floor) < 1e-9 else ""
        emit(f"  {head:24s}  {c['f_lo']:5.1f}-{c['f_hi']:5.1f} Hz  "
             f"p={c['p']:.4f}  p_fwer={c['p_fwer']:.4f}  "
             f"mass={c['mass']:+.1f}  ({direction}){floor}")
emit("-" * 68)


# ----------------------------------------------------------------------
# Per-contrast figures
# ----------------------------------------------------------------------
fmask = freqs <= FREQ_MAX_PLOT

def plot_contrast(ax, r):
    A, B = r["A"], r["B"]
    a = psd[:, pidx[A], :]
    b = psd[:, pidx[B], :]
    for ph, arr in [(A, a), (B, b)]:
        mn = arr.mean(0); se = arr.std(0, ddof=1) / np.sqrt(n_subj)
        ax.plot(freqs[fmask], mn[fmask], color=PHASE_COLORS[ph], lw=2, label=ph)
        ax.fill_between(freqs[fmask], (mn - se)[fmask], (mn + se)[fmask],
                        color=PHASE_COLORS[ph], alpha=0.22)
    ax.axhline(0, color="k", lw=0.6, alpha=0.5)
    # shade significant clusters
    ymin, ymax = ax.get_ylim()
    for c in r["clusters"]:
        if c["p"] < ALPHA_SIG:
            band = (freqs >= c["f_lo"]) & (freqs <= c["f_hi"])
            ax.fill_between(freqs, ymin, ymax, where=band & fmask,
                            color="gold", alpha=0.30, zorder=0)
            ax.text((c["f_lo"] + c["f_hi"]) / 2, ymax * 0.92,
                    f"p={c['p']:.3f}", ha="center", fontsize=8, fontweight="bold")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(freqs[fmask].min(), FREQ_MAX_PLOT)
    ax.set_title(f"{A} vs {B}", fontsize=11)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Power (dB re rest)")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.15)

# overview grid
n = len(results)
ncol = 2
nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(13, 3.4 * nrow))
for ax, r in zip(axes.ravel(), results):
    plot_contrast(ax, r)
for ax in axes.ravel()[n:]:
    ax.axis("off")
fig.suptitle("Subject-level CBPT (n=8) — STN-LFP phase spectra, dB re rest",
             fontsize=14, y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "CBPT_overview.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# individual PMBR + movement-vs-baseline figures at higher detail
for r in results:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plot_contrast(ax, r)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"CBPT_{r['A']}_vs_{r['B']}.png"),
                dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# Save machine-readable results
# ----------------------------------------------------------------------
json_out = {
    "n_subjects": n_subj,
    "subjects": subjects,
    "n_trials_per_subject": n_trials_subj.tolist(),
    "df": df,
    "cluster_threshold_t": float(thresh),
    "n_permutations": int(n_perm_exact),
    "p_floor": float(p_floor),
    "alpha_cluster": ALPHA_CLUSTER,
    "alpha_sig": ALPHA_SIG,
    "n_contrasts": n_contrasts,
    "contrasts": [],
}
for r in results:
    json_out["contrasts"].append({
        "A": r["A"], "B": r["B"], "family": r["family"],
        "clusters": [{k: c[k] for k in
                      ("f_lo", "f_hi", "n_bins", "mass", "t_peak", "p", "p_fwer", "sign")}
                     for c in r["clusters"]],
    })
with open(os.path.join(OUT_DIR, "clusters.json"), "w") as f:
    json.dump(json_out, f, indent=2)

with open(os.path.join(OUT_DIR, "summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"\nSaved figures + clusters.json + summary.txt to {OUT_DIR}")
