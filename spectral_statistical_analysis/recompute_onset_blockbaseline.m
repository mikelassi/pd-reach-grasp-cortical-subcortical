%% recompute_onset_blockbaseline.m
% =========================================================================
% Re-normalise the onset-locked LFP TF to the BLOCK pre-start rest baseline
% (the clean rest after the sync trigger, before any movement in the block),
% instead of the contaminated [-0.5, 0] s pre-movement window.
%
% Each trial is matched to its block's baseline_LFP (stored in the per-block
% *_LFP_trialsByRegionAndPhase.mat). The block baseline spectrum is the
% geometric-mean CWT power over the (edge-trimmed) rest segment, interpolated
% onto the epoch frequency axis. Trials are dB-normalised to that spectrum and
% median-averaged per subject.
%
% Output: onset_tf_blockbase.mat (-v7) — same fields as onset_tf.mat, ready
%         for CBPT_LFP_onset.py.
%
% Author: Michael Lassi
% =========================================================================
clear; clc;

EP_MAT = 'F:\Projects\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\Epochs\Epochs_allSubjects.mat';
BASE   = 'F:\Projects\Parkinson_ReachGrasp\Reprocessing';
SUB    = fullfile('Preprocessed','LFP');
OUT    = fullfile(fileparts(mfilename('fullpath')), 'onset_tf_blockbase.mat');

EDGE_TRIM_S = 0.2;   % trim each end of the baseline segment (CWT edge effects)

S = load(EP_MAT, 'EPOCHS'); EPOCHS = S.EPOCHS;
SUBJECTS = EPOCHS.params.subjects;
n_subj   = numel(SUBJECTS);

fs    = EPOCHS.(SUBJECTS{1}).fs;
n_on  = size(EPOCHS.(SUBJECTS{1}).onset.data, 1);
t_onset = EPOCHS.(SUBJECTS{1}).onset.times(:);

% Epoch filterbank (matches onset_tf.mat: 600 samples -> 66 freqs)
fb_ep = cwtfilterbank('SignalLength', n_on, 'SamplingFrequency', fs, ...
    'Wavelet','amor','FrequencyLimits',[1 80],'VoicesPerOctave',12);
fr = centerFrequencies(fb_ep);
n_fr = numel(fr);

per_subj_onset = nan(n_fr, n_on, n_subj);
base_quality   = nan(n_subj, 1);   % shortest block-baseline length (s)

for s = 1:n_subj
    subj = SUBJECTS{s};
    data = EPOCHS.(subj).onset.data;          % [n_on x n_trials]
    meta = EPOCHS.(subj).onset.meta;
    ch   = EPOCHS.(subj).ch_contra;
    blk  = [meta.block];
    n_tr = size(data, 2);

    p = fullfile(BASE, subj, SUB);
    mats = dir(fullfile(p, '*_LFP_trialsByRegionAndPhase.mat'));
    nblk = max(blk);

    % ---- per-block baseline spectrum on the epoch frequency axis ----
    base_spec = nan(n_fr, nblk);
    min_base_s = inf;
    for b = 1:nblk
        M  = load(fullfile(p, mats(b).name), 'baseline_LFP');
        bl = double(M.baseline_LFP(ch, :));
        tr = round(EDGE_TRIM_S * fs);
        if numel(bl) > 2*tr + 20, bl = bl(tr+1:end-tr); end
        min_base_s = min(min_base_s, numel(bl)/fs);
        fb_b = cwtfilterbank('SignalLength', numel(bl), 'SamplingFrequency', fs, ...
            'Wavelet','amor','FrequencyLimits',[1 80],'VoicesPerOctave',12);
        fr_b = centerFrequencies(fb_b);
        pw   = abs(wt(fb_b, bl)).^2;                      % [n_fr_b x len]
        spec = 10.^(mean(log10(max(pw, eps)), 2));        % geomean over time
        base_spec(:, b) = interp1(fr_b, spec, fr, 'linear', 'extrap');
    end
    base_quality(s) = min_base_s;

    % ---- per-trial CWT, dB re block baseline ----
    acc = nan(n_fr, n_on, n_tr);
    for t = 1:n_tr
        pw = abs(wt(fb_ep, data(:, t)')).^2;              % [n_fr x n_on]
        acc(:, :, t) = 10 * log10(max(pw, eps) ./ base_spec(:, blk(t)));
    end
    per_subj_onset(:, :, s) = median(acc, 3);
    fprintf('%-7s : %d trials, %d blocks, min block-baseline %.1f s\n', ...
        subj, n_tr, nblk, base_quality(s));
end

% COI mask (analytic Morlet e-folding) for the 600-sample epoch
coi_mask_onset = true(n_fr, n_on);
for fi = 1:n_fr
    cs = floor(fs * sqrt(2) / (2*pi*fr(fi)));
    coi_mask_onset(fi, 1:min(cs, n_on))            = false;
    coi_mask_onset(fi, max(1, n_on-cs+1):n_on)     = false;
end

subjects = SUBJECTS;
save(OUT, 'per_subj_onset', 't_onset', 'fr', 'subjects', ...
    'coi_mask_onset', 'base_quality', '-v7');
fprintf('\nSaved %s  [%d fr x %d t x %d subj]\n', OUT, n_fr, n_on, n_subj);
