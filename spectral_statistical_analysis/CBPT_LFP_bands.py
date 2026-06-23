"""
CBPT_LFP_bands.py
=========================================================================
A-priori CANONICAL-BAND analysis of STN-LFP phase spectra (subject-level,
exact sign-flip), with explicit multiple-comparison control.

Rationale
---------
The 1-80 Hz cluster test spends its multiple-comparison budget over 76
frequencies, so a narrow-band effect must out-compete the whole spectrum.
Restricting to the standard frequency bands -- which are specified a priori
from the STN/PD literature, NOT chosen from these data -- is a legitimate and
much more powerful test. To keep it honest, EVERY band x contrast test is
reported and corrected together (Benjamini-Hochberg FDR + Bonferroni); nothing
is cherry-picked.

Per cell: per-subject band-mean of the median (over trials) rest-normalised
spectrum -> exact subject-level sign-flip test (2^8 = 256 perms, floor 0.0039),
two-tailed. Direction and one-tailed p reported descriptively.

Secondary: individual beta-peak-aligned test for the movement-vs-baseline
contrasts. Each subject's beta peak frequency is taken from Rest_post (an
INDEPENDENT phase) to avoid circularity, then the movement ERD is measured in a
+/-2.5 Hz window around that subject-specific frequency.

Input : cbpt_input.mat
Output: CBPT_results/bands_heatmap.png, bands_table.csv, bands_summary.txt

Author: Michael Lassi
=========================================================================
"""

import os
import sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
IN_MAT = os.path.join(HERE, "cbpt_input.mat")
OUT_DIR = os.path.join(HERE, "CBPT_results")
os.makedirs(OUT_DIR, exist_ok=True)

BANDS = {
    "delta (1-4)": (1, 4), "theta (4-8)": (4, 8), "alpha (8-13)": (8, 13),
    "low-beta (13-20)": (13, 20), "high-beta (20-30)": (20, 30), "gamma (30-60)": (30, 60),
}
CONTRASTS = [
    ("Rest_pre", "Reach"), ("Rest_pre", "Grasp"), ("Rest_pre", "Pull"),
    ("Rest_pre", "Rest_post"),
    ("Reach", "Grasp"), ("Grasp", "Pull"), ("Reach", "Pull"),
]
ALPHA = 0.05

# ----------------------------------------------------------------------
def _cellstr(a):
    return [str(x.ravel()[0]) if isinstance(x, np.ndarray) and x.size else str(x)
            for x in np.asarray(a).ravel()]

m = sio.loadmat(IN_MAT)
freqs = m["freqs"].ravel().astype(float)
phase_names = _cellstr(m["phase_names"])
psd = np.asarray(m["psd_subj_phase"], dtype=float)   # [S,P,F] medians
pi = {p: i for i, p in enumerate(phase_names)}
S = psd.shape[0]

# exhaustive sign matrix
signs = np.array(np.meshgrid(*([[1, -1]] * S))).T.reshape(-1, S).astype(float)
n_perm = signs.shape[0]
p_floor = 1.0 / n_perm


def signflip(x, tail=0):
    """Exact sign-flip test of mean(x)=0. Returns (t_obs, p)."""
    obs = x.mean() / (x.std(ddof=1) / np.sqrt(S) + 1e-12)
    perm = (signs * x[None]).mean(1) / ((signs * x[None]).std(1, ddof=1) / np.sqrt(S) + 1e-12)
    if tail == 0:
        p = np.mean(np.abs(perm) >= abs(obs) - 1e-9)
    elif tail < 0:
        p = np.mean(perm <= obs + 1e-9)
    else:
        p = np.mean(perm >= obs - 1e-9)
    return obs, p


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        val = p[i] * n / (rank + 1)
        prev = min(prev, val)
        q[i] = prev
    return np.minimum(q, 1.0)


# ----------------------------------------------------------------------
# Pooled "Movement" phase = mean of Reach/Grasp/Pull (one test, more power)
psd_move = psd[:, [pi["Reach"], pi["Grasp"], pi["Pull"]], :].mean(1)   # [S,F]

def band_diff(specA, specB, lo, hi):
    fm = (freqs >= lo) & (freqs <= hi)
    return specA[:, fm].mean(1) - specB[:, fm].mean(1)

# ----- PRIMARY confirmatory family: pre-specified, DIRECTIONAL -----
# Directions fixed a priori from the STN/PD reach-to-grasp literature.
# One-tailed -> sign-flip floor = 1/256 = 0.0039, so Bonferroni across the
# small primary family can still reach 0.05 at n=8.
PRIMARY = [
    # well-established directions -> one-tailed; beta split into low/high
    ("PMBR low-beta rebound",   psd[:, pi["Rest_post"], :],
     psd[:, pi["Rest_pre"], :], "low-beta (13-20)",  +1),
    ("PMBR high-beta rebound",  psd[:, pi["Rest_post"], :],
     psd[:, pi["Rest_pre"], :], "high-beta (20-30)", +1),
    ("Movement low-beta ERD",   psd_move,
     psd[:, pi["Rest_pre"], :], "low-beta (13-20)",  -1),
    ("Movement high-beta ERD",  psd_move,
     psd[:, pi["Rest_pre"], :], "high-beta (20-30)", -1),
    # direction not firmly established in STN -> two-tailed (honest)
    ("Movement theta modulation", psd_move,
     psd[:, pi["Rest_pre"], :], "theta (4-8)", 0),
]
# Exploratory/secondary (NOT pre-specified; emerged from the data) — reported
# separately, never counted in the confirmatory family.
SECONDARY = [
    ("Pull>Reach high-beta build-up", psd[:, pi["Pull"], :],
     psd[:, pi["Reach"], :], "high-beta (20-30)", +1),
]
PRIMARY_BANDS = {"low-beta (13-20)": (13, 20), "high-beta (20-30)": (20, 30),
                 "theta (4-8)": (4, 8)}

band_names = list(BANDS.keys())
nb, nc = len(band_names), len(CONTRASTS)
T = np.zeros((nc, nb)); P = np.zeros((nc, nb)); D = np.zeros((nc, nb)); P1 = np.zeros((nc, nb))

for ci, (A, B) in enumerate(CONTRASTS):
    for bi, bn in enumerate(band_names):
        lo, hi = BANDS[bn]
        fm = (freqs >= lo) & (freqs <= hi)
        d = psd[:, pi[A], :][:, fm].mean(1) - psd[:, pi[B], :][:, fm].mean(1)
        t, p2 = signflip(d, 0)
        _, p1 = signflip(d, -1 if d.mean() < 0 else 1)
        T[ci, bi], P[ci, bi], D[ci, bi], P1[ci, bi] = t, p2, d.mean(), p1

q = bh_fdr(P.ravel()).reshape(P.shape)
bonf = np.minimum(P * P.size, 1.0)

# ----------------------------------------------------------------------
rep = []
def emit(s=""):
    print(s); rep.append(s)

emit("#" * 78)
emit("PRIMARY (confirmatory) — pre-specified directional hypotheses, one-tailed")
emit(f"   sign-flip floor (1-tailed) = {1.0/n_perm:.4f};  Bonferroni m = {len(PRIMARY)}")
emit("#" * 78)
prim_p = []
prim_rows = []
for (name, specA, specB, bname, tail) in PRIMARY:
    lo, hi = PRIMARY_BANDS[bname]
    d = band_diff(specA, specB, lo, hi)
    t, p1 = signflip(d, tail)
    prim_p.append(p1)
    prim_rows.append((name, bname, tail, d.mean(), t, p1))
prim_bonf = np.minimum(np.array(prim_p) * len(PRIMARY), 1.0)
for (name, bname, tail, dm, t, p1), pb in zip(prim_rows, prim_bonf):
    verdict = "CONFIRMED" if pb < ALPHA else ("trend" if p1 < ALPHA else "n.s.")
    arrow = "(1-tail A>B)" if tail > 0 else ("(1-tail A<B)" if tail < 0 else "(2-tailed)")
    emit(f"  {name:34s} {bname:12s} {arrow}")
    emit(f"      d={dm:+5.2f}dB  t={t:+5.2f}  p={p1:.4f}  p_Bonf={pb:.4f}  -> {verdict}")
emit("")
emit("Secondary (exploratory, NOT in confirmatory family — needs confirmation):")
for (name, specA, specB, bname, tail) in SECONDARY:
    lo, hi = PRIMARY_BANDS[bname]
    d = band_diff(specA, specB, lo, hi)
    t, p1 = signflip(d, tail)
    emit(f"  {name:34s} {bname:12s}  d={d.mean():+5.2f}dB t={t:+5.2f} p={p1:.4f} (uncorrected)")
emit("")

emit("=" * 78)
emit("EXPLORATORY — full canonical-band grid (two-tailed, hypothesis-generating)")
emit(f"Family = {nb} bands x {nc} contrasts = {P.size} tests; BH-FDR reported.")
emit("NOTE: at n=8 the 2-tailed sign-flip floor is %.4f, so BH-FDR across %d tests"
     % (2.0 / n_perm, P.size))
emit("      cannot reach q<0.05 for any cell — this is why pre-specification (above)")
emit("      is required, not optional. The grid below is exploratory only.")
emit("=" * 78)
hdr = f"{'contrast':22s}" + "".join(f"{b.split(' ')[0]:>9s}" for b in band_names)
emit(hdr)
for ci, (A, B) in enumerate(CONTRASTS):
    line = f"{A+' vs '+B:22s}"
    for bi in range(nb):
        mark = "+" if q[ci, bi] < ALPHA else ("." if P[ci, bi] < ALPHA else " ")
        line += f"{T[ci,bi]:>7.1f}{mark} "
    emit(line)
emit("")
emit("Legend: '+' FDR q<0.05 (survives correction) | '.' raw p<0.05 only")
emit("")
emit("Significant after FDR correction:")
for ci, (A, B) in enumerate(CONTRASTS):
    for bi in range(nb):
        if q[ci, bi] < ALPHA:
            dirn = f"{A}>{B}" if D[ci, bi] > 0 else f"{A}<{B}"
            emit(f"  {A:9s} vs {B:9s} | {band_names[bi]:13s} "
                 f"t={T[ci,bi]:+5.2f} p={P[ci,bi]:.4f} q={q[ci,bi]:.4f} "
                 f"d={D[ci,bi]:+.2f}dB ({dirn})")

# ---- secondary: individual beta-peak-aligned movement ERD ----
emit("")
emit("-" * 78)
emit("Secondary: individual beta-peak-aligned (peak from Rest_post; +/-2.5 Hz)")
bm = (freqs >= 13) & (freqs <= 30); fb = freqs[bm]
peakf = np.array([fb[np.argmax(psd[s, pi["Rest_post"], bm])] for s in range(S)])
emit("  per-subject beta peak (Hz): " + ", ".join(f"{p:.1f}" for p in peakf))
for (A, B) in [("Rest_pre", "Reach"), ("Rest_pre", "Grasp"),
               ("Rest_pre", "Pull"), ("Reach", "Pull")]:
    d = np.array([psd[s, pi[A], (freqs >= peakf[s]-2.5) & (freqs <= peakf[s]+2.5)].mean()
                  - psd[s, pi[B], (freqs >= peakf[s]-2.5) & (freqs <= peakf[s]+2.5)].mean()
                  for s in range(S)])
    t, p2 = signflip(d, 0); _, p1 = signflip(d, -1 if d.mean() < 0 else 1)
    emit(f"  {A:9s} vs {B:9s}: d={d.mean():+5.2f}dB t={t:+5.2f} p2t={p2:.4f} p1t={p1:.4f}")
emit("-" * 78)

# ----------------------------------------------------------------------
# Heatmap of t-values with significance overlay
fig, ax = plt.subplots(figsize=(9, 5))
vlim = np.nanmax(np.abs(T))
im = ax.imshow(T, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
ax.set_xticks(range(nb)); ax.set_xticklabels(band_names, rotation=25, ha="right")
ax.set_yticks(range(nc)); ax.set_yticklabels([f"{A} vs {B}" for A, B in CONTRASTS])
for ci in range(nc):
    for bi in range(nb):
        mark = "**" if q[ci, bi] < ALPHA else ("*" if P[ci, bi] < ALPHA else "")
        ax.text(bi, ci, f"{T[ci,bi]:.1f}{mark}", ha="center", va="center",
                fontsize=8, color="black")
cb = fig.colorbar(im, ax=ax); cb.set_label("sign-flip t  (A vs B)")
ax.set_title("Canonical-band STN-LFP contrasts (n=8)\n** FDR q<0.05   * raw p<0.05",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "bands_heatmap.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# CSV
import csv
with open(os.path.join(OUT_DIR, "bands_table.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["contrast", "band", "t", "p_2tail", "q_FDR", "p_Bonf", "mean_dB", "p_1tail"])
    for ci, (A, B) in enumerate(CONTRASTS):
        for bi, bn in enumerate(band_names):
            w.writerow([f"{A} vs {B}", bn, f"{T[ci,bi]:.3f}", f"{P[ci,bi]:.4f}",
                        f"{q[ci,bi]:.4f}", f"{bonf[ci,bi]:.4f}",
                        f"{D[ci,bi]:.3f}", f"{P1[ci,bi]:.4f}"])

with open(os.path.join(OUT_DIR, "bands_summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(rep))

print(f"\nSaved bands_heatmap.png + bands_table.csv + bands_summary.txt to {OUT_DIR}")
