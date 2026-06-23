"""
CBPT_LFP_onset.py
=========================================================================
Event-locked (NON-warped) time-frequency cluster test around movement onset.

Motivation
----------
The movement beta-ERD is a transient locked to movement onset. The warped
phase-average smears that brief dip across the whole Reach window, which is
why it reads as a weak ~-0.6 dB effect. Here we keep real time (epoch
[-0.5, +1.0] s around onset A, dB re the [-0.5, 0] s pre-movement baseline)
and let a short, deep ERD form a cluster in time instead of being averaged out.

Design (identical inference machinery to the warped TF-plane test)
  X = per-subject onset-locked TF, dB re pre-move rest  [n_subj, n_freq, n_time]
  One-sample cluster test of X != 0; exhaustive subject-level sign-flip
  (2^8 = 256, floor 0.0039); 2-D (freq x time) lattice adjacency.
  Threshold: fixed t_{7} (two-tailed) and TFCE.
  Frequency restricted 4-80 Hz a priori; cone-of-influence pixels masked out.

Input : onset_tf.mat (export_onset_tf.m)
Output: CBPT_results/onset_fixed.png, onset_tfce.png, onset_clusters.json

Author: Michael Lassi
=========================================================================
"""

import os, sys, json
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
IN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "onset_tf.mat")
if not os.path.isabs(IN):
    IN = os.path.join(HERE, IN)
PREFIX = sys.argv[2] if len(sys.argv) > 2 else "onset"
OUT_DIR = os.path.join(HERE, "CBPT_results")
os.makedirs(OUT_DIR, exist_ok=True)

FREQ_LO, FREQ_HI = 4, 80
ALPHA = 0.05
CAXIS = 3.0

m = sio.loadmat(IN)
fr = m["fr"].ravel().astype(float)                 # [n_fr] descending
t = m["t_onset"].ravel().astype(float)             # [n_time] seconds
P = np.asarray(m["per_subj_onset"], float)         # [n_fr, n_time, n_subj]
coi = np.asarray(m["coi_mask_onset"], bool)        # [n_fr, n_time] true = valid

X = np.transpose(P, (2, 0, 1))                     # [S, n_fr, n_time]
n_subj, n_fr, n_time = X.shape

fmask = (fr >= FREQ_LO) & (fr <= FREQ_HI)
fr_u = fr[fmask]
Xb = X[:, fmask, :].copy()                         # [S, nf, nt]
coi_b = coi[fmask, :]                              # [nf, nt]
nf = fmask.sum()

# zero-out cone-of-influence pixels (conservative: no effect there -> cannot
# seed or extend a cluster). Baseline region is ~0 dB anyway.
Xb[:, ~coi_b] = 0.0

df = n_subj - 1
thr = st.t.ppf(1 - ALPHA / 2, df)
n_perm = 2 ** n_subj

print("=" * 70, flush=True)
print("Event-locked (onset) one-sample TF CBPT (vs pre-movement baseline)", flush=True)
print("=" * 70, flush=True)
print(f"X = [{n_subj} subj, {nf} freq ({FREQ_LO}-{FREQ_HI} Hz), {n_time} time "
      f"({t[0]:+.2f}..{t[-1]:+.2f} s)]", flush=True)
print(f"Exhaustive sign-flip: {n_perm} perms (floor {1/n_perm:.4f}); "
      f"fixed |t|>{thr:.3f} and TFCE", flush=True)
print("=" * 70, flush=True)


def run(threshold):
    T_obs, clusters, pv, _ = permutation_cluster_1samp_test(
        Xb, threshold=threshold, n_permutations=n_perm, tail=0,
        adjacency=None, out_type="mask", seed=0, verbose=False)
    sig = np.zeros((nf, n_time), bool); info = []
    for cl, p in zip(clusters, pv):
        if p < ALPHA:
            sig |= cl
            fi, ti = np.where(cl)
            mass = T_obs[cl].sum()
            info.append(dict(
                f_lo=float(fr_u[fi].min()), f_hi=float(fr_u[fi].max()),
                t_lo=float(t[ti.min()]), t_hi=float(t[ti.max()]),
                n_pix=int(cl.sum()), mass=float(mass), p=float(p),
                sign="ERS(+)" if mass > 0 else "ERD(-)"))
    return T_obs, sig, info


def plot(sig, info, tag, fname):
    grand = Xb.mean(0)
    order = np.argsort(fr_u); fr_s = fr_u[order]
    g = grand[order]; s = sig[order]
    fig, ax = plt.subplots(figsize=(9, 5))
    pcm = ax.pcolormesh(t, fr_s, g, cmap="RdBu_r", vmin=-CAXIS, vmax=CAXIS, shading="auto")
    ax.contour(t, fr_s, s.astype(float), levels=[0.5], colors="k", linewidths=1.6)
    ax.axvline(0, color="k", ls="--", lw=1.5)
    ax.text(0, fr_s[-1]*0.98, "onset", rotation=90, va="top", ha="right", fontsize=8)
    ax.axhline(13, color="0.4", ls=":", lw=0.8); ax.axhline(30, color="0.4", ls=":", lw=0.8)
    ax.set_yscale("log"); ax.set_yticks([4, 8, 13, 20, 30, 50, 80])
    ax.set_yticklabels([4, 8, 13, 20, 30, 50, 80])
    ax.set_xlabel("Time from movement onset (s)"); ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"STN-LFP onset-locked CBPT ({tag}) — outline = cluster p<{ALPHA}\n"
                 f"n={n_subj}, {len(info)} significant cluster(s)")
    cb = fig.colorbar(pcm, ax=ax); cb.set_label("dB re pre-move rest")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=200, bbox_inches="tight")
    plt.close(fig)


results = {}
for threshold, tag, fname in [(thr, "fixed-threshold", f"{PREFIX}_fixed.png"),
                              (dict(start=0, step=0.2), "TFCE", f"{PREFIX}_tfce.png")]:
    T_obs, sig, info = run(threshold)
    plot(sig, info, tag, fname)
    results[tag] = info
    print(f"\n## {tag}: {len(info)} significant cluster(s)", flush=True)
    for c in sorted(info, key=lambda d: d["p"]):
        print(f"  {c['sign']}  {c['f_lo']:4.1f}-{c['f_hi']:4.1f} Hz  "
              f"t {c['t_lo']:+.2f}..{c['t_hi']:+.2f} s  p={c['p']:.4f}  "
              f"mass={c['mass']:+.0f}  ({c['n_pix']} px)", flush=True)

with open(os.path.join(OUT_DIR, f"{PREFIX}_clusters.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved {PREFIX}_fixed.png, {PREFIX}_tfce.png, {PREFIX}_clusters.json to {OUT_DIR}", flush=True)
