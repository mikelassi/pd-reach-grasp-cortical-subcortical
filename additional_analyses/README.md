# Cortex – STN – Muscle interaction (multivariate spectral connectivity)

Quantifies how the **motor cortex**, the **subthalamic nucleus (STN)** and the
**contralateral muscles** interact during a reach-to-grasp task in 8 STN-DBS
Parkinson's patients (OFF medication), using EEG (cortex), bipolar STN LFP and
surface EMG recorded synchronously at 400 Hz.

All three connectivity measures are derived from a **single multitaper
cross-spectral density (CSD) matrix** per condition, so they are mutually
consistent:

| Measure | Question it answers | Directed? | Conditioned? |
|---|---|---|---|
| **Coherence** | Are two nodes linearly coupled at frequency *f*? | no | no |
| **Partial coherence** | Is the coupling *direct*, or routed via other nodes? | no | yes (on the rest of the network) |
| **Granger causality** (non-parametric, Wilson factorisation) | Which node *drives* which? | yes | pairwise |

## Nodes

* **Cortical ROIs** (groups of scalp EEG channels, ROI signal = mean of
  per-channel z-scored traces), laterality-corrected so *contralateral* always
  means the hemisphere opposite the moving hand:
  `M1_contra`, `PM_contra` (premotor), `S1par_contra` (post-central/parietal),
  `SMA` (mesial midline), `M1_ipsi` (control).
* **STN** — contralateral bipolar LFP channel (`ch_contra` from `Epoching.m`).
* **Muscles** — contralateral **IOD, Triceps, Deltoid**, band-passed 20–180 Hz,
  50/100 Hz notch, and **full-wave rectified** (standard CMC input; demodulates
  motor-unit firing so the common ~15–30 Hz descending drive appears in the EMG
  spectrum — Mima & Hallett 1999).

Laterality (from `preprocessing/BatchEMGProcessing_v2.m`): left-hand movers
(`wue02/03/10/11`) → contralateral = right hemisphere; right-hand movers
(`wue05/06/07/09`) → contralateral = left hemisphere.

### IMPORTANT — EMG event alignment (fixed)

The kinematic events (A–F) are detected on the **wrist velocity inside the
EMG_KIN file** (`A02_events_segmentation.m`) and stored in *that file's* sample
coordinates in `02_Kinematics/Events/*_kinematic_block.mat`. The Preprocessed
EEG/LFP were **start-cropped** during preprocessing, so their embedded event
latencies live in a *different, cropped* timeline (LFP files are 3–18 s shorter
than the EMG_KIN files). Injecting LFP-timeline event latencies into the
uncropped EMG file — as the lab's `CMC_analysis.m` does — misaligns the EMG by
**several seconds**, which destroys any genuine corticomuscular coherence.
This pipeline therefore epochs **EMG with its own native events** (from the
`*_kinematic_block.mat`) and EEG/LFP with their own embedded events; epochs are
then matched across modalities **by trial number** (`_T<n>`), not by position.

## Conditions (windows reused from the lab CMC / TF pipeline)

* **Movement** `[0, +1.0] s` around movement onset (event A)
* **Pull** `[0, +0.875] s` around pull onset (event D)
* **Baseline** `[-0.5, 0] s` pre-movement rest (event A)

## Files

| File | Purpose |
|---|---|
| `cortex_stn_muscle_connectivity.m` | **Main pipeline.** Loads EEG/STN/EMG, epochs on the same trials, builds the node set, computes the full multivariate network + focused `M1–STN–Muscle` triads per condition, statistics, saves results, draws figures. |
| `cs_csd_multitaper.m` | Multitaper (Slepian) cross-spectral density matrix. |
| `cs_coherence.m` | Magnitude-squared coherence from a CSD. |
| `cs_partial_coherence.m` | Partial coherence from the inverse CSD (precision matrix). |
| `cs_wilson.m` | Wilson minimum-phase spectral factorisation (engine for non-parametric Granger). |
| `cs_granger.m` | Pairwise non-parametric spectral Granger causality (Geweke). |
| `cs_bandavg.m` | Band-averaging helper. |
| `make_figures.m` | Renders all figures from the saved `.mat` (decoupled from computation). |
| `validate_connectivity.m` | **Validation harness** — proves every estimator on synthetic ground truth. Run this first. |
| `check_results.m` | Prints a group-level numeric summary of the saved results. |
| `diag_cmc_investigate.m` | Group-level controls: file-length mismatch, EMG movement-locking, intermuscular coherence, full-spectrum CMC. |
| `diag_hold.m` | CMC across movement phases (reach A / grasp-hold C / pull D). |
| `diag_laplacian.m` | CMC with ROI-mean vs single-electrode vs surface (Hjorth) Laplacian referencing. |
| `inspect_data.m`, `diag_*.m` | Other scratch diagnostics from development (data layout, EMG rectification, ROI choice, EMG label selection, native-event alignment). |

## How to run

```matlab
cd additional_analyses
validate_connectivity            % must print "ALL TESTS PASSED"
cortex_stn_muscle_connectivity   % full analysis (set QUICK_TEST=true for a 1-subject smoke test)
make_figures                     % (optional) re-draw figures from saved results
check_results                    % (optional) numeric summary
```

Data path is auto-detected (`F:\Projects\Parkinson_ReachGrasp\Reprocessing`,
falling back to the legacy `H:\...`). Requires EEGLAB (for `.set` loading) and
the Signal Processing Toolbox (`dpss`, `butter`, `filtfilt`).

### Key options (top of the main script)
* `QUICK_TEST` — single subject, no surrogates (pipeline smoke test).
* `N_SURR` — trial-shuffle surrogates for triad-Granger significance (0 to skip).
* `EMG_RECTIFY` — full-wave rectify EMG before connectivity (default true).
* `NW`, `K` — multitaper time-bandwidth and taper count.

## Outputs → `RESULTS_final/CortexSTNMuscle/`

* `CortexSTNMuscle_results.mat` — all per-subject and grand-average results.
* `network_matrices_beta.png` — beta-band coherence / partial coherence /
  Granger matrices over the whole network, per condition.
* `triad_coh_spectra.png`, `triad_pcoh_spectra.png` — triad coherence and
  partial-coherence spectra (M1–STN, STN–Muscle, M1–Muscle), per condition.
* `triad_granger_spectra.png` — directed Granger spectra (dashed = surrogate
  95% threshold when `N_SURR>0`).
* `crosssystem_beta_summary.png` — headline bar chart of beta-band coherence for
  all cross-system links vs the 95% confidence limit.

## Validation

`validate_connectivity.m` checks each estimator against analytic ground truth:

1. **Wilson factorisation** reconstructs an exact MVAR spectrum to ~1e-11 and
   recovers the true noise covariance to ~1e-16.
2. **Coherence** recovers a known analytic value (mixed common-source model).
3. **Partial coherence** collapses an *indirect* link in an X→Y→Z chain
   (coh(X,Z)≈0.48 → partial(X,Z|Y)≈0.00) while keeping the direct links.
4. **Granger** recovers the correct *direction* in a unidirectional AR system
   (GC ratio > 500:1).

## Statistical inference

* **Coherence / partial coherence:** analytic 95% confidence limit
  `1 - alpha^(1/(K·n_trials - 1))` (multitaper effective DOF; Halliday 1995),
  drawn on the spectra and summary bar chart.
* **Granger causality:** trial-shuffle surrogates (`N_SURR`) break cross-node
  trial pairing to build an empirical null; the 95th percentile is overlaid as a
  dashed threshold on the directed-Granger figure.

## Findings (grand average, n = 8)

* A robust **cortico-subthalamic beta interaction** (~20–30 Hz):
  `M1_contra–STN`, `SMA–STN` and `PM_contra–STN` coherence clearly exceed the
  95% limit, peaking in high beta; **directed Granger is bidirectional**
  (`M1→STN ≈ STN→M1`), consistent with the cortico-basal-ganglia beta loop in PD.
* **Corticomuscular and subthalamo-muscular coherence are weak** (peaks ~0.02–0.03,
  marginally above the single-bin 95% limit at best) in this OFF-medication,
  dynamic reach-to-grasp paradigm.
* The **intermuscular** common drive is, by contrast, **strong** (low-frequency
  coherence ≈ 0.35) — a positive control proving the EMG is clean and carries a
  shared descending drive; the weak *cortico*-muscular coupling is therefore a
  genuine property of the task/state, not a signal-quality problem.

### Investigation of the weak corticomuscular coherence

This was investigated exhaustively (scripts `diag_*.m`):

1. **EMG alignment** — the dominant artifact. Before the fix (injected
   cropped-LFP events) the EMG was misaligned by 3–18 s and CMC was pure noise;
   the native-event fix roughly doubled–tripled low-frequency cortico-muscular
   coherence and made the EMG envelope movement-locked.
2. **Frequency** — checked 2–100 Hz (not only beta): cortico-muscular coupling is
   largest at *low frequency* (common/movement drive) and in *high beta*
   (STN–muscle, Marsden 2001), not classic 13–30 Hz beta.
3. **Movement phase** — tested reach (A), grasp-hold (C) and pull (D); CMC is
   weak in all (`diag_hold.m`).
4. **EEG referencing** — ROI mean vs single electrode vs surface (Hjorth)
   Laplacian; none recovers strong CMC (`diag_laplacian.m`).

**Conclusion (literature-validated).** Strong corticomuscular coherence is *not*
expected here: Hirschmann et al. (2013, NeuroImage) show M1–muscle coherence is
*strongly reduced during repetitive/dynamic movement* (vs static contraction) in
PD; CMC is primarily an isometric-hold phenomenon and is further reduced in PD
OFF medication. Marsden, Brown et al. (2001, Brain) report STN–EMG coherence that
is high-beta (28–44 Hz) and *transient* (time-locked bursts), which a static
1-s window averages out. The strong, robust interaction in this dataset is the
**cortico-subthalamic beta loop**, with muscle coupling weak — the expected
signature for OFF-medication PD performing a phasic reach-to-grasp.

> Caution: partial coherence *among the cortical ROIs themselves* is affected by
> volume conduction / spatial collinearity and should not be over-interpreted.
> The scientific focus is the non-collinear cross-system links (cortex–STN,
> cortex–muscle, STN–muscle), where partial coherence is reliable.

## References
Wilson (1972); Geweke (1982); Halliday et al. (1995); Mima & Hallett (1999);
Marsden, Brown et al. (2001, *Brain*); Dahlhaus (2000); Dhamala, Rangarajan &
Ding (2008); Boonstra & Breakspear (2012); Hirschmann et al. (2013, *NeuroImage*).
