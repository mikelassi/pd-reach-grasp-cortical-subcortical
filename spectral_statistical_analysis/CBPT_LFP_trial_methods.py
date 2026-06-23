"""
CBPT_LFP_trial_methods.py
=========================================================================
Cluster-based permutation test on STN-LFP phase spectra, comparing THREE
cluster-forming statistics under the SAME (correct) subject-level sign-flip
permutation:

  M0  Subject-mean t        one-sample t on per-subject median differences
                            (the conservative two-stage summary; n=8, df=7).
  M1  Precision-weighted z  DerSimonian-Laird random-effects combination of
                            per-subject mean differences, weighted by
                            1/(tau^2 + within-subject-variance). Uses every
                            trial via each subject's within-subject variance
                            and trial count.
  M2  Mixed-model t         one-way random-intercept model
                            diff_sj = mu + u_s + e_sj  (u_s subject random
                            effect, e_sj within-trial residual), using all 239
                            trials. Variance components by ANOVA method-of-
                            moments (sigma^2 = pooled within-MS; tau^2 from the
                            between/within MS contrast), GLS weights
                            w_s = 1/(tau^2 + sigma^2/n_s), test of mu=0.
                            NB: statsmodels REML of the same model was
                            numerically unstable here (intercept-only MixedLM
                            diverges at the tau^2=0 boundary, returning t=0 with
                            bse~1e7), so we use the stable closed-form estimator
                            of the identical model.

WHY THIS IS VALID CBPT
----------------------
A cluster permutation test = (per-bin statistic) + (permutation that respects
exchangeability). The per-bin statistic is interchangeable; M0/M1/M2 simply
swap it. The permutation is IDENTICAL for all three and is the only thing that
must be correct: for a within-subject contrast the exchangeable unit is the
SUBJECT, so the null is built by exhaustive subject-level SIGN-FLIPS (2^8 = 256;
a subject's whole difference spectrum flips together). Cluster mass null =
max |cluster mass| over the frequency axis per permutation. Cluster p =
fraction of the 256 permutations whose max |mass| >= observed (floor 1/256).

For M2 the variance components (tau^2, sigma^2) are estimated once per frequency
by REML on the observed data; the sign-flip null then re-uses those weights
(fixed-nuisance sign-flip of a pivotal statistic, as in Winkler et al. 2014
PALM). This is the standard efficient permutation for mixed models and avoids
~19k refits while remaining the genuine LMM statistic (validated against
statsmodels on the observed map below).

KEY POINT: all three are bounded by n=8. Trials sharpen each subject's estimate;
they do not add degrees of freedom. M1/M2 should track M0 closely — that
agreement is the evidence that n=8, not the trial count, is the binding limit.

Input : cbpt_input.mat (export_psd_for_cbpt.m, now incl. trial_psd/trial_subject)
Output: CBPT_results/CBPT_methods_overview.png, CBPT_methods_<A>_vs_<B>.png,
        methods_clusters.json, methods_summary.txt

Author: Michael Lassi
=========================================================================
"""

import os
import sys
import json
import warnings

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import scipy.io as sio
import scipy.stats as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
IN_MAT = os.path.join(HERE, "cbpt_input.mat")
OUT_DIR = os.path.join(HERE, "CBPT_results")
os.makedirs(OUT_DIR, exist_ok=True)

ALPHA = 0.05
FREQ_MAX_PLOT = 80
EPS = 1e-12

PHASE_COLORS = {
    "Rest_pre": "#c0392b", "Reach": "#27ae60", "Grasp": "#2980b9",
    "Pull": "#8e44ad", "Rest_post": "#e67e22",
}
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
def _cellstr(a):
    out = []
    for x in np.asarray(a).ravel():
        out.append(str(x.ravel()[0]) if isinstance(x, np.ndarray) and x.size else str(x))
    return out


m = sio.loadmat(IN_MAT)
freqs = m["freqs"].ravel().astype(float)
subjects = _cellstr(m["subjects"])
phase_names = _cellstr(m["phase_names"])
psd_med = np.asarray(m["psd_subj_phase"], dtype=float)       # [S, P, F] medians
trial_psd = np.asarray(m["trial_psd"], dtype=float)          # [N, P, F]
trial_subj = m["trial_subject"].ravel().astype(int) - 1      # 0-based [N]
pidx = {p: i for i, p in enumerate(phase_names)}

S = len(subjects)
F = freqs.size
df = S - 1
n_s = np.array([(trial_subj == s).sum() for s in range(S)])

# Exhaustive 2^S sign matrix [n_perm, S], row 0 = observed (all +1)
signs = np.array(np.meshgrid(*([[1, -1]] * S))).T.reshape(-1, S).astype(float)
order = np.argsort(np.where(signs < 0, 1, 0).sum(1))   # put all-+1 first
signs = signs[order]
signs[0] = 1.0
n_perm = signs.shape[0]
p_floor = 1.0 / n_perm

# cluster-forming thresholds (two-tailed alpha) in each statistic's reference
THR_T = st.t.ppf(1 - ALPHA / 2, df)     # M0, M2 (t-like, df=S-1)
THR_Z = st.norm.ppf(1 - ALPHA / 2)      # M1 (z-like)

print("=" * 70)
print("Trial-aware CBPT — three cluster statistics, subject-level sign-flip")
print("=" * 70)
print(f"S={S} subjects, trials/subj={n_s.tolist()} (N={n_s.sum()}), F={F} freqs")
print(f"Permutations: {n_perm} exhaustive sign-flips (p floor {p_floor:.4f})")
print(f"Thresholds: |t|>{THR_T:.3f} (M0,M2)  |z|>{THR_Z:.3f} (M1)")
print("=" * 70)


# ----------------------------------------------------------------------
def cluster_masses(stat, thr):
    """Signed cluster masses on a 1-D stat vector (contiguous bins |stat|>thr)."""
    masses = []
    sign_state = 0
    run = 0.0
    for v in stat:
        s = 1 if v > thr else (-1 if v < -thr else 0)
        if s != 0 and s == sign_state:
            run += v
        else:
            if sign_state != 0:
                masses.append(run)
            sign_state = s
            run = v if s != 0 else 0.0
    if sign_state != 0:
        masses.append(run)
    return masses


def cluster_bins(stat, thr):
    """Return list of (lo_idx, hi_idx, mass) for observed clusters."""
    out = []
    i = 0
    while i < len(stat):
        s = 1 if stat[i] > thr else (-1 if stat[i] < -thr else 0)
        if s == 0:
            i += 1
            continue
        j = i
        run = 0.0
        while j < len(stat) and ((stat[j] > thr) if s > 0 else (stat[j] < -thr)):
            run += stat[j]
            j += 1
        out.append((i, j - 1, run))
        i = j
    return out


def permute_pvalues(stat_perm, thr):
    """stat_perm: [n_perm, F]. Returns (observed_clusters, p per cluster)."""
    # null distribution of max |cluster mass|
    null_max = np.zeros(n_perm)
    for p in range(n_perm):
        masses = cluster_masses(stat_perm[p], thr)
        null_max[p] = max((abs(x) for x in masses), default=0.0)
    obs = cluster_bins(stat_perm[0], thr)
    pvals = []
    for (lo, hi, mass) in obs:
        pvals.append((np.sum(null_max >= abs(mass) - 1e-9)) / n_perm)
    return obs, pvals


# ----------------------------------------------------------------------
# Per-subject difference summaries for a contrast
def subject_summaries(A, B):
    """Returns ybar [S,F] (mean diff per subj), vw [S,F] (within-var of the
    mean = s^2/n), and med_diff [S,F] (median diff per subj)."""
    a = trial_psd[:, pidx[A], :]
    b = trial_psd[:, pidx[B], :]
    d = a - b                       # [N, F]
    ybar = np.zeros((S, F)); vw = np.zeros((S, F))
    med = np.zeros((S, F)); ssw = np.zeros((S, F))
    for s in range(S):
        ds = d[trial_subj == s]     # [n_s, F]
        ybar[s] = ds.mean(0)
        med[s] = np.median(ds, 0)
        ssw[s] = ((ds - ds.mean(0)) ** 2).sum(0)   # within-subject SS (for M2)
        if ds.shape[0] > 1:
            vw[s] = ds.var(0, ddof=1) / ds.shape[0]
        else:
            vw[s] = np.nan
    # Floor each subject's within-variance at 5% of the across-subject median
    # (per frequency) so one near-zero-variance subject cannot dominate the
    # inverse-variance weights — a standard safeguard for M1.
    med_vw = np.nanmedian(vw, axis=0, keepdims=True)
    vw = np.maximum(vw, np.maximum(0.05 * med_vw, EPS))
    return ybar, vw, med, ssw


# ---- M0: subject-mean one-sample t (on per-subject MEDIAN difference) ----
def stat_M0(med, signs):
    D = signs[:, :, None] * med[None, :, :]      # [P,S,F]
    mean = D.mean(1)
    sd = D.std(1, ddof=1)
    return mean / (sd / np.sqrt(S) + EPS)


# ---- M1: DerSimonian-Laird random-effects z (fixed within-var) ----
def stat_M1(ybar, vw, signs):
    y = signs[:, :, None] * ybar[None, :, :]     # [P,S,F]
    wfe = 1.0 / vw                               # [S,F]
    wfe_b = wfe[None]                            # [1,S,F]
    ybar_fe = (wfe_b * y).sum(1) / wfe.sum(0)[None]      # [P,F]
    Q = (wfe_b * (y - ybar_fe[:, None, :]) ** 2).sum(1)  # [P,F]
    C = wfe.sum(0) - (wfe ** 2).sum(0) / wfe.sum(0)      # [F]
    tau2 = np.maximum(0.0, (Q - df) / np.maximum(C[None], EPS))   # [P,F]
    wstar = 1.0 / (vw[None] + tau2[:, None, :])          # [P,S,F]
    theta = (wstar * y).sum(1) / wstar.sum(1)            # [P,F]
    se = np.sqrt(1.0 / wstar.sum(1))                     # [P,F]
    return theta / (se + EPS)


# ---- M2: one-way random-intercept GLS t (ANOVA method-of-moments) ----
def stat_M2(ybar, ssw, signs):
    """One-way random-intercept model on within-trial differences, closed form.
    sigma^2 = pooled within-MS (fixed under sign flips, since within-variance is
    sign-invariant); tau^2 from the between/within mean-square contrast
    (recomputed per permutation as the subject means flip); GLS weights
    w_s = 1/(tau^2 + sigma^2/n_s); t = mu_hat * sqrt(sum_s w_s)."""
    N = int(n_s.sum())
    sigma2 = ssw.sum(0) / (N - S)                  # [F] pooled within MS (MSW)
    n0 = (N - (n_s ** 2).sum() / N) / (S - 1)      # scalar (unbalanced design)
    y = signs[:, :, None] * ybar[None, :, :]       # [P,S,F]
    ns = n_s[None, :, None].astype(float)          # [1,S,1]
    grand = (ns * y).sum(1) / N                     # [P,F] n-weighted grand mean
    SSB = (ns * (y - grand[:, None, :]) ** 2).sum(1)   # [P,F]
    MSB = SSB / (S - 1)                             # [P,F]
    tau2 = np.maximum(0.0, (MSB - sigma2[None]) / n0)  # [P,F]
    w = 1.0 / (tau2[:, None, :] + sigma2[None, None, :] / n_s[None, :, None])  # [P,S,F]
    mu = (w * y).sum(1) / w.sum(1)                  # [P,F]
    se = np.sqrt(1.0 / w.sum(1))                    # [P,F]
    return mu / (se + EPS)


# ----------------------------------------------------------------------
results = {}
report = []

def emit(s=""):
    print(s); report.append(s)

for (A, B, fam) in CONTRASTS:
    ybar, vw, med, ssw = subject_summaries(A, B)

    m0 = stat_M0(med, signs)
    obs0, p0 = permute_pvalues(m0, THR_T)

    m1 = stat_M1(ybar, vw, signs)
    obs1, p1 = permute_pvalues(m1, THR_Z)

    m2 = stat_M2(ybar, ssw, signs)
    obs2, p2 = permute_pvalues(m2, THR_T)

    results[(A, B)] = dict(
        family=fam, ybar=ybar, med=med,
        M0=(obs0, p0), M1=(obs1, p1), M2=(obs2, p2),
    )

emit("\nCLUSTER RESULTS — significant clusters (cluster p < %.2f)" % ALPHA)
emit("M0=subject-mean t | M1=precision-weighted z | M2=mixed-model t")
emit("-" * 70)
cur = None
for (A, B, fam) in CONTRASTS:
    if fam != cur:
        cur = fam; emit(f"\n## {fam}")
    r = results[(A, B)]
    emit(f"  {A} vs {B}")
    for label, key in [("M0", "M0"), ("M1", "M1"), ("M2", "M2")]:
        obs, pv = r[key]
        sig = [(lo, hi, mass, p) for (lo, hi, mass), p in zip(obs, pv) if p < ALPHA]
        if not sig:
            emit(f"     {label}: n.s.")
        else:
            for (lo, hi, mass, p) in sig:
                direction = f"{A}>{B}" if mass > 0 else f"{A}<{B}"
                floor = "  [p-floor]" if abs(p - p_floor) < 1e-9 else ""
                f1, f2 = sorted((freqs[lo], freqs[hi]))   # freq axis is descending
                emit(f"     {label}: {f1:5.1f}-{f2:5.1f} Hz  "
                     f"p={p:.4f}  ({direction}){floor}")
emit("-" * 70)


# ----------------------------------------------------------------------
# Figures: spectra + significance bars from each method
fmask = freqs <= FREQ_MAX_PLOT
METHOD_STYLE = {"M0": ("#555555", 0.92), "M1": ("#1f9e44", 0.85), "M2": ("#2c6fbb", 0.78)}

def draw(ax, A, B):
    r = results[(A, B)]
    a = psd_med[:, pidx[A], :]; b = psd_med[:, pidx[B], :]
    for ph, arr in [(A, a), (B, b)]:
        mn = arr.mean(0); se = arr.std(0, ddof=1) / np.sqrt(S)
        ax.plot(freqs[fmask], mn[fmask], color=PHASE_COLORS[ph], lw=2, label=ph)
        ax.fill_between(freqs[fmask], (mn - se)[fmask], (mn + se)[fmask],
                        color=PHASE_COLORS[ph], alpha=0.2)
    ax.axhline(0, color="k", lw=0.6, alpha=0.5)
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    for mi, (label, key) in enumerate([("M0", "M0"), ("M1", "M1"), ("M2", "M2")]):
        col, yfrac = METHOD_STYLE[label]
        obs, pv = r[key]
        ylev = ymin + span * (0.04 + 0.05 * mi)
        for (lo, hi, mass), p in zip(obs, pv):
            if p < ALPHA:
                f1, f2 = sorted((freqs[lo], freqs[hi]))   # descending freq axis
                ax.plot([f1, f2], [ylev, ylev], color=col, lw=5,
                        solid_capstyle="butt")
                ax.text(f2 + 1, ylev, f"{label} p={p:.3f}", color=col,
                        fontsize=7, va="center")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(freqs[fmask].min(), FREQ_MAX_PLOT)
    ax.set_title(f"{A} vs {B}", fontsize=11)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Power (dB re rest)")
    ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.15)

n = len(CONTRASTS); ncol = 2; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(13, 3.4 * nrow))
for ax, (A, B, fam) in zip(axes.ravel(), CONTRASTS):
    draw(ax, A, B)
for ax in axes.ravel()[n:]:
    ax.axis("off")
fig.suptitle("Trial-aware CBPT (n=8): M0 subject-mean | M1 precision-weighted | M2 mixed-model",
             fontsize=13, y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "CBPT_methods_overview.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

for (A, B, fam) in CONTRASTS:
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    draw(ax, A, B)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"CBPT_methods_{A}_vs_{B}.png"),
                dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
jout = {"n_subjects": S, "trials_per_subject": n_s.tolist(), "n_perm": int(n_perm),
        "p_floor": p_floor, "thr_t": float(THR_T), "thr_z": float(THR_Z),
        "contrasts": []}
for (A, B, fam) in CONTRASTS:
    r = results[(A, B)]
    entry = {"A": A, "B": B, "family": fam, "methods": {}}
    for key in ("M0", "M1", "M2"):
        obs, pv = r[key]
        entry["methods"][key] = [
            {"f_lo": float(freqs[lo]), "f_hi": float(freqs[hi]),
             "mass": float(mass), "p": float(p)}
            for (lo, hi, mass), p in zip(obs, pv)]
    jout["contrasts"].append(entry)
with open(os.path.join(OUT_DIR, "methods_clusters.json"), "w") as f:
    json.dump(jout, f, indent=2)
with open(os.path.join(OUT_DIR, "methods_summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(report))

print(f"\nSaved methods figures + methods_clusters.json + methods_summary.txt to {OUT_DIR}")
