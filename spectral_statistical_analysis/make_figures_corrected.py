"""
make_figures_corrected.py
=========================================================================
Polished, paneled summary figure of the CORRECTED LFP spectral statistics:
  A  5-phase rest-normalised PSD (mean +/- SEM across subjects)
  B  PMBR contrast (Rest_pre vs Rest_post) with the subject-level CBPT
     significant cluster shaded
  C  Reach-phase high-beta ERD spectrum (vs block rest), beta band marked
  D  Reach high-beta (20-30 Hz) ERD per subject (beta-peak group): the
     significant group result (p = 0.016, 7/7)

Data:
  cbpt_input.mat  (psd_subj_phase: dB re pre-movement rest, 5 phases)
  reach_psd.mat   (reach/grip/rest PSDs, dB re BLOCK rest)

Author: Michael Lassi
=========================================================================
"""
import os, sys
import numpy as np
import scipy.io as sio
import scipy.stats as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "CBPT_results"); os.makedirs(OUT, exist_ok=True)

# ---- aesthetics ----
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 240, "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.labelsize": 11, "legend.fontsize": 9, "legend.frameon": False,
    "axes.grid": True, "grid.alpha": 0.18, "grid.linewidth": 0.6,
})
PH = ["Rest_pre", "Reach", "Grasp", "Pull", "Rest_post"]
COL = {"Rest_pre": "#7f8c8d", "Reach": "#27ae60", "Grasp": "#2980b9",
       "Pull": "#8e44ad", "Rest_post": "#e67e22"}
SIG = "#f1c40f"        # significant-cluster gold
HIBETA = "#e74c3c"

def _cell(a): return [str(x.ravel()[0]) if isinstance(x, np.ndarray) and x.size else str(x)
                      for x in np.asarray(a).ravel()]

# ---- load ----
m = sio.loadmat(os.path.join(HERE, "cbpt_input.mat"))
fr = m["freqs"].ravel().astype(float)
phn = _cell(m["phase_names"]); pidx = {p: i for i, p in enumerate(phn)}
psd = np.asarray(m["psd_subj_phase"], float)          # [S,5,F] dB re pre-rest
S = psd.shape[0]

r = sio.loadmat(os.path.join(HERE, "reach_psd.mat"))
frR = r["freqs"].ravel().astype(float)
R = np.asarray(r["reach_psd"], float); B = np.asarray(r["rest_psd"], float)
HASPK = r["HASPK"].ravel().astype(bool); subs = _cell(r["SUBJECTS"])

# ---- subject-level sign-flip cluster for PMBR (Rest_pre vs Rest_post) ----
signs = np.array(np.meshgrid(*([[1, -1]] * S))).T.reshape(-1, S).astype(float)
def cluster_freqs(A, Bp, lo=4, hi=45):
    fm = (fr >= lo) & (fr <= hi); f = fr[fm]
    X = psd[:, pidx[A], fm] - psd[:, pidx[Bp], fm]
    df = S - 1; thr = st.t.ppf(0.975, df)
    def tmap(Y): return Y.mean(0) / (Y.std(0, ddof=1) / np.sqrt(S) + 1e-12)
    def cl(tv):
        out = []; i = 0
        while i < len(tv):
            s = 1 if tv[i] > thr else (-1 if tv[i] < -thr else 0)
            if s == 0: i += 1; continue
            j = i; mass = 0.0
            while j < len(tv) and ((tv[j] > thr) if s > 0 else (tv[j] < -thr)): mass += tv[j]; j += 1
            out.append((i, j - 1, mass)); i = j
        return out
    obs = cl(tmap(X)); null = np.zeros(len(signs))
    for k in range(len(signs)):
        c = cl(tmap(signs[k][:, None] * X)); null[k] = max((abs(x[2]) for x in c), default=0)
    res = []
    for (a, b, mass) in obs:
        p = np.mean(null >= abs(mass) - 1e-9)
        res.append((f[a], f[b], p))
    return f, res
fpm, pmbr_cl = cluster_freqs("Rest_pre", "Rest_post")
pmbr_sig = [(a, b, p) for (a, b, p) in pmbr_cl if p < 0.05]

# ---- band helper (sign-flip) ----
def band_test(X, lo, hi, Sset, tail=0):
    fm = (frR >= lo) & (frR <= hi); x = X[Sset][:, fm].mean(1); n = len(Sset)
    sg = np.array(np.meshgrid(*([[1, -1]] * n))).T.reshape(-1, n).astype(float)
    o = x.mean() / (x.std(ddof=1) / np.sqrt(n) + 1e-12)
    pr = (sg * x[None]).mean(1) / ((sg * x[None]).std(1, ddof=1) / np.sqrt(n) + 1e-12)
    p = (np.mean(np.abs(pr) >= abs(o) - 1e-9) if tail == 0 else np.mean(pr <= o + 1e-9))
    return x, x.mean(), o, p

# =====================================================================
fig = plt.figure(figsize=(13.5, 9))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.24,
                       left=0.07, right=0.97, top=0.93, bottom=0.08)

def panel_label(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="top")

# ---- Panel A: 5-phase PSD overview ----
axA = fig.add_subplot(gs[0, 0])
fmA = (fr >= 2) & (fr <= 45); fA = fr[fmA]
for ph in PH:
    mn = psd[:, pidx[ph], fmA].mean(0); se = psd[:, pidx[ph], fmA].std(0, ddof=1) / np.sqrt(S)
    axA.fill_between(fA, mn - se, mn + se, color=COL[ph], alpha=0.18)
    axA.plot(fA, mn, color=COL[ph], lw=2.2, label=ph.replace("_", " "))
axA.axhline(0, color="k", lw=0.7)
for (a, b, p) in pmbr_sig:
    axA.axvspan(a, b, color=SIG, alpha=0.20, zorder=0)
axA.set_xlim(2, 45); axA.set_xlabel("Frequency (Hz)"); axA.set_ylabel("Power (dB re pre-move rest)")
axA.set_title("Kinematic-phase spectra (grand average ± SEM)")
axA.legend(ncol=2, loc="upper right")
panel_label(axA, "A")

# ---- Panel B: PMBR contrast ----
axB = fig.add_subplot(gs[0, 1])
for ph in ["Rest_pre", "Rest_post"]:
    mn = psd[:, pidx[ph], fmA].mean(0); se = psd[:, pidx[ph], fmA].std(0, ddof=1) / np.sqrt(S)
    axB.fill_between(fA, mn - se, mn + se, color=COL[ph], alpha=0.22)
    axB.plot(fA, mn, color=COL[ph], lw=2.4, label=ph.replace("_", " "))
axB.axhline(0, color="k", lw=0.7)
for (a, b, p) in pmbr_sig:
    axB.axvspan(a, b, color=SIG, alpha=0.28, zorder=0)
    axB.text((a + b) / 2, axB.get_ylim()[1] * 0.86, f"cluster\np = {p:.3f}",
             ha="center", fontsize=9, fontweight="bold")
axB.set_xlim(2, 45); axB.set_xlabel("Frequency (Hz)"); axB.set_ylabel("Power (dB re pre-move rest)")
axB.set_title("Post-movement beta rebound (subject-level CBPT)")
axB.legend(loc="upper right")
panel_label(axB, "B")

# ---- Panel C: Reach high-beta ERD spectrum (vs block rest) ----
axC = fig.add_subplot(gs[1, 0])
keep = [i for i in range(S) if HASPK[i]]
X = 10 * np.log10(R / B)                         # reach vs block rest, dB
fmC = (frR >= 4) & (frR <= 45); fC = frR[fmC]
mn = X[keep][:, fmC].mean(0); se = X[keep][:, fmC].std(0, ddof=1) / np.sqrt(len(keep))
axC.axvspan(20, 30, color=HIBETA, alpha=0.10, zorder=0, label="high-β (20–30 Hz)")
axC.fill_between(fC, mn - se, mn + se, color=COL["Reach"], alpha=0.22)
axC.plot(fC, mn, color=COL["Reach"], lw=2.4)
axC.axhline(0, color="k", lw=0.7)
axC.set_xlim(4, 45); axC.set_xlabel("Frequency (Hz)")
axC.set_ylabel("Reach power (dB re block rest)")
axC.set_title("Reach-phase desynchronisation (n=7, beta-peak)")
axC.legend(loc="lower right")
panel_label(axC, "C")

# ---- Panel D: Reach high-beta ERD per subject (the significant result) ----
axD = fig.add_subplot(gs[1, 1])
vals, md, t, p = band_test(X, 20, 30, keep, tail=0)
bp = axD.boxplot(vals, vert=True, widths=0.45, patch_artist=True, showfliers=False,
                 positions=[1], medianprops=dict(color="k", lw=2),
                 boxprops=dict(facecolor=HIBETA, alpha=0.22, edgecolor="k"),
                 whiskerprops=dict(color="k"), capprops=dict(color="k"))
rng = np.random.default_rng(0)
jit = 1 + (rng.random(len(vals)) - 0.5) * 0.22
axD.scatter(jit, vals, s=70, color=HIBETA, edgecolor="k", zorder=3, linewidth=0.6)
axD.scatter([1], [vals.mean()], marker="D", s=90, color="k", zorder=4, label="mean")
axD.axhline(0, color="k", lw=0.9, ls="--")
axD.set_xticks([1]); axD.set_xticklabels(["Reach"])
axD.set_ylabel("High-β power (dB re block rest)")
axD.set_title("High-β reach ERD — group result")
ymax = max(vals.max(), 0.3)
axD.text(1.32, vals.mean(), f"mean = {md:+.2f} dB\np = {p:.3f}\n{(vals<0).sum()}/{len(vals)} ERD",
         va="center", fontsize=10, fontweight="bold")
axD.set_xlim(0.5, 2.0)
panel_label(axD, "D")

fig.suptitle("STN-LFP reach-to-grasp — corrected subject-level spectral statistics (n=8; n=7 with beta-peak criterion)",
             fontsize=13.5, fontweight="bold", y=0.985)
out = os.path.join(OUT, "LFP_corrected_summary.png")
fig.savefig(out, bbox_inches="tight")
print("PMBR significant clusters:", [(round(a,1),round(b,1),round(p,4)) for a,b,p in pmbr_sig])
print(f"Reach high-beta ERD (n={len(keep)}): mean={md:+.2f} dB, t={t:+.2f}, p={p:.4f}")
print("Saved", out)
