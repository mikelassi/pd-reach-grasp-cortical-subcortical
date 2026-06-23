"""
CBPT_LFP_TFplane.py
=========================================================================
2-D time-frequency-plane cluster-based permutation test on the warped STN-LFP
TF maps (Phase 2).

Why this is the right test for the beta ERD
-------------------------------------------
The band/marginal tests collapse each phase to one number, which dilutes a
time-localised event-related desynchronisation. Here we keep the full warped
(frequency x movement-time) plane and ask, with proper multiple-comparison
control, WHERE power departs from the pre-movement rest baseline (0 dB). A
brief, strong beta drop around movement onset can now form a cluster in time
instead of being averaged away.

Design
------
  X = per-subject warped TF map, dB re rest        [n_subj, n_freq, n_time]
  One-sample cluster test of X != 0 (rest baseline = 0 dB by construction).
  Permutation : exhaustive subject-level sign-flip (2^8 = 256; floor 0.0039).
  Adjacency   : 2-D lattice over (frequency, warped-time) — adjacency=None.
  Threshold   : (a) fixed cluster-forming t_{7} (two-tailed), and
                (b) TFCE (threshold-free), which is more sensitive to broad,
                    low-amplitude effects like a distributed ERD.
  Frequency restricted to 4-80 Hz a priori (below ~4 Hz the Morlet support
  exceeds the short warped phases, so those estimates are unreliable).

Input : TF_SUBJ_allSubjects.mat  (from TF_analysis_LFP_v2.m, -v7)
Output: CBPT_results/TFplane_fixed.png, TFplane_tfce.png, TFplane_clusters.json

Author: Michael Lassi
=========================================================================
"""

import os
import sys
import json
import numpy as np
import scipy.io as sio
import scipy.stats as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mne.stats import permutation_cluster_1samp_test

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
TF_MAT = r"F:\Projects\Parkinson_ReachGrasp\Reprocessing\RESULTS_compare\v2\Baseline_01_median\TF_SUBJ_allSubjects.mat"
OUT_DIR = os.path.join(HERE, "CBPT_results")
os.makedirs(OUT_DIR, exist_ok=True)

FREQ_LO, FREQ_HI = 4, 80
ALPHA = 0.05
CAXIS = 3.0   # dB display limit

# ----------------------------------------------------------------------
m = sio.loadmat(TF_MAT, struct_as_record=False, squeeze_me=True)
TF = m["TF_SUBJ"]
fr = np.asarray(TF.fr, dtype=float).ravel()                 # [n_fr] descending
tpct = np.asarray(TF.time_axis_pct, dtype=float).ravel()    # [n_time]
tf_db = np.asarray(TF.tf_db, dtype=float)                    # [n_fr, n_time, n_subj]
phase_pct = np.asarray(TF.phase_pct, dtype=float).ravel()   # [onset, grasp, pull, offset]

# -> [n_subj, n_fr, n_time]
X = np.transpose(tf_db, (2, 0, 1))
n_subj, n_fr, n_time = X.shape

# restrict frequency band a priori
fmask = (fr >= FREQ_LO) & (fr <= FREQ_HI)
fr_u = fr[fmask]
Xb = X[:, fmask, :]                                          # [S, nf, nt]
nf = fmask.sum()

df = n_subj - 1
thr = st.t.ppf(1 - ALPHA / 2, df)
n_perm = 2 ** n_subj

print("=" * 70)
print("TF-plane one-sample CBPT (vs rest baseline = 0 dB)")
print("=" * 70)
print(f"X = [{n_subj} subj, {nf} freq ({FREQ_LO}-{FREQ_HI} Hz), {n_time} warped-time]")
print(f"Exhaustive sign-flip: {n_perm} perms (floor {1/n_perm:.4f})")
print(f"Fixed cluster-forming threshold |t|>{thr:.3f}; also TFCE")
print("=" * 70)


def run(threshold, tag):
    T_obs, clusters, cluster_pv, _ = permutation_cluster_1samp_test(
        Xb, threshold=threshold, n_permutations=n_perm, tail=0,
        adjacency=None, out_type="mask", seed=0, verbose=False)
    sig = np.zeros((nf, n_time), dtype=bool)
    info = []
    for cl, pv in zip(clusters, cluster_pv):
        if pv < ALPHA:
            sig |= cl
            fidx, tidx = np.where(cl)
            fr_in = fr_u[fidx]
            mass = T_obs[cl].sum()
            info.append(dict(
                f_lo=float(fr_in.min()), f_hi=float(fr_in.max()),
                t_lo=float(tpct[tidx.min()]), t_hi=float(tpct[tidx.max()]),
                n_pix=int(cl.sum()), mass=float(mass), p=float(pv),
                sign="ERS(+)" if mass > 0 else "ERD(-)"))
    return T_obs, sig, info


def plot(T_obs, sig, info, tag, fname):
    grand = Xb.mean(0)                            # [nf, nt]
    # sort ascending in freq for display
    order = np.argsort(fr_u)
    fr_s = fr_u[order]
    g = grand[order]; s = sig[order]
    fig, ax = plt.subplots(figsize=(11, 5))
    pcm = ax.pcolormesh(tpct, fr_s, g, cmap="RdBu_r", vmin=-CAXIS, vmax=CAXIS,
                        shading="auto")
    # outline significant clusters
    ax.contour(tpct, fr_s, s.astype(float), levels=[0.5], colors="k", linewidths=1.6)
    for x in phase_pct:
        ax.axvline(x, color="0.25", ls="--", lw=1)
    for x, lab in zip(phase_pct, ["onset", "grasp", "pull", "offset"]):
        ax.text(x, fr_s[-1] * 0.98, lab, rotation=90, va="top", ha="right",
                fontsize=8, color="0.25")
    ax.set_yscale("log")
    ax.set_yticks([4, 8, 13, 20, 30, 50, 80])
    ax.set_yticklabels([4, 8, 13, 20, 30, 50, 80])
    ax.set_xlabel("Warped movement time (%)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"STN-LFP TF-plane CBPT ({tag}) — black outline = cluster p<{ALPHA}\n"
                 f"n={n_subj}, {len(info)} significant cluster(s)")
    cb = fig.colorbar(pcm, ax=ax); cb.set_label("dB re rest")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=200, bbox_inches="tight")
    plt.close(fig)


results = {}
for threshold, tag, fname in [
    (thr, "fixed-threshold", "TFplane_fixed.png"),
    (dict(start=0, step=0.2), "TFCE", "TFplane_tfce.png"),
]:
    T_obs, sig, info = run(threshold, tag)
    plot(T_obs, sig, info, tag, fname)
    results[tag] = info
    print(f"\n## {tag}: {len(info)} significant cluster(s)")
    for c in sorted(info, key=lambda d: d["p"]):
        print(f"  {c['sign']}  {c['f_lo']:4.1f}-{c['f_hi']:4.1f} Hz  "
              f"time {c['t_lo']:+5.0f}..{c['t_hi']:+5.0f}%  "
              f"p={c['p']:.4f}  mass={c['mass']:+.0f}  ({c['n_pix']} px)")

with open(os.path.join(OUT_DIR, "TFplane_clusters.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved TFplane_fixed.png, TFplane_tfce.png, TFplane_clusters.json to {OUT_DIR}")
