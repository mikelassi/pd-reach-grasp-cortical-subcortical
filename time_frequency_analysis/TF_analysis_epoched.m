%% TF_analysis_epoched.m
% =========================================================================
% Time-Frequency analysis of fixed-length STN LFP epochs using CWT.
%
% Input:
%   Epochs_allSubjects.mat  (from Epoching.m)
%   Contralateral-hemisphere channel, fixed-length epochs:
%     MovOnset : [-0.5, +1.0] s around event A   (600 samples at 400 Hz)
%     Pull     : [-0.3, +0.875] s around event D  (470 samples at 400 Hz)
%
% Key differences from TF_analysis_LFP_v2.m:
%   1. No time-warping: epochs are already fixed length, direct CWT per trial.
%   2. Single shared baseline: the pre-movement REST window [-0.5, 0] s
%      from the onset epoch is used to normalise BOTH epoch types, giving a
%      consistent rest-referenced dB scale across conditions.
%   3. COI masking: the cone-of-influence is computed and masked on all TF
%      plots — critical for short epochs where low frequencies are affected
%      at the epoch edges.
%   4. Smoothing applied for display only (does not affect saved matrices).
%
% Normalisation:
%   dB = 10 * log10(power / baseline)
%   baseline = geometric mean (mean in log domain) of power in [-0.5, 0] s
%              from onset epochs → median across trials → per subject.
%
% Outputs (saved to RESULTS_DIR/TF_epoched/):
%   TF_grandAvg.png / .fig
%   TF_perSubject_MovOnset.png / .fig
%   TF_perSubject_Pull.png / .fig
%   TF_allSubjects.mat
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
MAT_FILE    = 'H:\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\Epochs\Epochs_allSubjects.mat';
RESULTS_DIR = 'H:\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\TF_epoched';

% CWT — same settings as TF_analysis_LFP_v2 for comparability
FREQ_LIMITS     = [1, 80];
VOICES_PER_OCT  = 12;
WAVELET_TYPE    = 'amor';   % analytic Morlet

% Display
FREQ_PLOT_MIN   = 4;    % Hz — below ~4 Hz COI is significant for short epochs
FREQ_PLOT_MAX   = 80;
CAXIS_DB        = 5;    % ± dB for normalised colormap
TF_SMOOTH_SIGMA = [1, 6];   % [freq, time] samples — display only

% Baseline window (s) from onset epoch (rest before movement)
BASELINE_WIN    = [-0.5, 0];

%% ===== INIT =====
fprintf('Loading epochs...\n');
load(MAT_FILE, 'EPOCHS');
SUBJECTS = EPOCHS.params.subjects;
n_subj   = numel(SUBJECTS);

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

% Colormap: blue–white–red
n_c = 256;
ramp = linspace(0,1,n_c/2)';
bwr  = [[ramp, ramp, ones(n_c/2,1)]; [ones(n_c/2,1), flip(ramp), flip(ramp)]];

% Smoothing lambda
if all(TF_SMOOTH_SIGMA > 0)
    tf_smooth = @(M) imgaussfilt(double(M), TF_SMOOTH_SIGMA);
else
    tf_smooth = @(M) double(M);
end

% Frequency band reference lines for plots
BANDS = struct('name',{'θ','α','β','lγ'}, ...
               'freq',{8, 13, 30, 60}, ...
               'col', {[.6 .4 0],[.4 .0 .8],[.0 .5 .0],[.8 .0 .0]});

%% ===== PRE-PASS: build single shared CWT filterbank =====
% Pull epochs (470 samp) are zero-padded to onset length (600 samp) so
% both epoch types share an identical frequency axis and baseline reference.
% The padded region falls in the COI for low frequencies and is masked out.
fs_ref  = EPOCHS.(SUBJECTS{1}).fs;
n_onset = size(EPOCHS.(SUBJECTS{1}).onset.data, 1);   % 600
n_pull  = size(EPOCHS.(SUBJECTS{1}).pull.data,  1);   % 470

fb = cwtfilterbank('SignalLength', n_onset, 'SamplingFrequency', fs_ref, ...
    'Wavelet', WAVELET_TYPE, 'FrequencyLimits', FREQ_LIMITS, 'VoicesPerOctave', VOICES_PER_OCT);

fr    = centerFrequencies(fb);   % [n_fr × 1], descending (shared by both conditions)
n_fr  = numel(fr);

t_onset = EPOCHS.(SUBJECTS{1}).onset.times;
t_pull  = EPOCHS.(SUBJECTS{1}).pull.times;

% COI masks — analytic Morlet e-folding time: sqrt(2)/(2*pi*f) seconds.
% Mask is computed for the valid (non-padded) duration of each epoch.
coi_mask_onset = true(n_fr, n_onset);
coi_mask_pull  = true(n_fr, n_pull);
for fi = 1:n_fr
    coi_samp = floor(fs_ref * sqrt(2) / (2 * pi * fr(fi)));
    coi_mask_onset(fi, 1:min(coi_samp,n_onset))           = false;
    coi_mask_onset(fi, max(1,n_onset-coi_samp+1):n_onset) = false;
    coi_mask_pull(fi,  1:min(coi_samp,n_pull))            = false;
    coi_mask_pull(fi,  max(1,n_pull-coi_samp+1):n_pull)   = false;
end

fprintf('CWT: %d log-spaced frequencies [%.1f–%.1f Hz], %d VPO\n', ...
    n_fr, fr(end), fr(1), VOICES_PER_OCT);

%% ===== SUBJECT LOOP =====
all_tf_onset  = zeros(n_fr, n_onset, n_subj);   % grand-avg accumulator
all_tf_pull   = zeros(n_fr, n_pull,  n_subj);
all_bl        = zeros(n_fr, n_subj);             % per-subject baseline

for s = 1:n_subj
    subj = SUBJECTS{s};
    if ~isfield(EPOCHS, subj)
        warning('Subject %s not found in EPOCHS — skipping.', subj);
        continue;
    end

    fs        = EPOCHS.(subj).fs;
    data_on   = EPOCHS.(subj).onset.data;   % [600 × n_trials]
    data_pu   = EPOCHS.(subj).pull.data;    % [470 × n_trials]
    n_t_on    = EPOCHS.(subj).onset.n;
    n_t_pu    = EPOCHS.(subj).pull.n;

    fprintf('\n===== Subject: %s  (onset=%d, pull=%d trials) =====\n', ...
        subj, n_t_on, n_t_pu);

    % ---------------------------------------------------------------
    % STEP 1: CWT on onset epochs → extract per-trial rest baseline
    % ---------------------------------------------------------------
    bl_mask = t_onset >= BASELINE_WIN(1) & t_onset < BASELINE_WIN(2);

    pow_onset_all = zeros(n_fr, n_onset, n_t_on);
    trial_baselines = zeros(n_fr, n_t_on);

    for t = 1:n_t_on
        cwt_t = wt(fb, data_on(:,t)');               % [n_fr × n_onset]
        pow_t = abs(cwt_t).^2;
        pow_onset_all(:,:,t) = pow_t;

        % Geometric-mean baseline (mean in log domain — avoids Jensen bias)
        trial_baselines(:,t) = 10.^(mean(log10(max(pow_t(:, bl_mask), eps)), 2));
    end

    % Subject baseline: median across trials (robust to burst artifacts)
    subj_baseline = median(trial_baselines, 2);   % [n_fr × 1]
    all_bl(:,s)   = subj_baseline;

    % ---------------------------------------------------------------
    % STEP 2: dB normalise onset epochs using subject baseline
    % ---------------------------------------------------------------
    pow_onset_db = 10 * log10(max(pow_onset_all, eps) ./ subj_baseline);
    tf_onset_subj = median(pow_onset_db, 3);   % [n_fr × n_onset]

    % ---------------------------------------------------------------
    % STEP 3: CWT on pull epochs, normalise with SAME subject baseline
    % ---------------------------------------------------------------
    pow_pull_all = zeros(n_fr, n_pull, n_t_pu);

    for t = 1:n_t_pu
        % Zero-pad pull epoch to n_onset so the shared filterbank can be used
        sig_padded = [data_pu(:,t)', zeros(1, n_onset - n_pull)];
        cwt_t = wt(fb, sig_padded);                  % [n_fr × n_onset]
        pow_pull_all(:,:,t) = abs(cwt_t(:, 1:n_pull)).^2;   % trim padding
    end

    pow_pull_db   = 10 * log10(max(pow_pull_all, eps) ./ subj_baseline);
    tf_pull_subj  = median(pow_pull_db, 3);    % [n_fr × n_pull]

    % ---------------------------------------------------------------
    % STEP 4: Store for grand average
    % ---------------------------------------------------------------
    all_tf_onset(:,:,s) = tf_onset_subj;
    all_tf_pull(:,:,s)  = tf_pull_subj;

    % ---------------------------------------------------------------
    % Per-subject TF stored for saving
    % ---------------------------------------------------------------
    EPOCHS.(subj).tf_onset = tf_onset_subj;
    EPOCHS.(subj).tf_pull  = tf_pull_subj;
    EPOCHS.(subj).baseline = subj_baseline;

    fprintf('  Done.\n');
end

%% ===== GRAND AVERAGE =====
grand_tf_onset = mean(all_tf_onset, 3);   % mean across subjects (one median per subject)
grand_tf_pull  = mean(all_tf_pull,  3);
sem_tf_onset   = std(all_tf_onset,  0, 3) / sqrt(n_subj);
sem_tf_pull    = std(all_tf_pull,   0, 3) / sqrt(n_subj);

%% ===== FIGURE 1: Grand Average TF =====
figure('Name','Grand Average TF','Position',[50 50 1300 530],'Visible','off');

ax1 = subplot(1,2,1);
plot_tf_panel(ax1, t_onset, fr, grand_tf_onset, coi_mask_onset, bwr, ...
    CAXIS_DB, FREQ_PLOT_MIN, FREQ_PLOT_MAX, tf_smooth, BANDS, ...
    sprintf('Grand Average — Movement Onset  (n=%d)', n_subj), 'A (onset)', 0);

ax2 = subplot(1,2,2);
plot_tf_panel(ax2, t_pull, fr, grand_tf_pull, coi_mask_pull, bwr, ...
    CAXIS_DB, FREQ_PLOT_MIN, FREQ_PLOT_MAX, tf_smooth, BANDS, ...
    sprintf('Grand Average — Pull Onset  (n=%d)', n_subj), 'D (pull)', 0);

sgtitle('STN LFP Time-Frequency Power (dB re pre-movement rest)', 'FontSize', 14);
set(gcf,'Units','pixels','Position',[50 50 1300 530]);
print(fullfile(RESULTS_DIR,'TF_grandAvg'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'TF_grandAvg.fig'));

%% ===== FIGURE 2: Per-Subject TF — Movement Onset =====
figure('Name','Per-Subject TF — MovOnset','Position',[50 50 1600 800],'Visible','off');
for s = 1:n_subj
    subj = SUBJECTS{s};
    ax = subplot(2,4,s);
    if ~isfield(EPOCHS, subj) || ~isfield(EPOCHS.(subj), 'tf_onset')
        title(ax, [subj ' — missing']); continue
    end
    plot_tf_panel(ax, t_onset, fr, EPOCHS.(subj).tf_onset, coi_mask_onset, bwr, ...
        CAXIS_DB, FREQ_PLOT_MIN, FREQ_PLOT_MAX, tf_smooth, BANDS, subj, 'A', 0);
end
sgtitle('Per-Subject TF — Movement Onset (dB re rest)', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1600 800]);
print(fullfile(RESULTS_DIR,'TF_perSubject_MovOnset'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'TF_perSubject_MovOnset.fig'));

%% ===== FIGURE 3: Per-Subject TF — Pull =====
figure('Name','Per-Subject TF — Pull','Position',[50 50 1600 800],'Visible','off');
for s = 1:n_subj
    subj = SUBJECTS{s};
    ax = subplot(2,4,s);
    if ~isfield(EPOCHS, subj) || ~isfield(EPOCHS.(subj), 'tf_pull')
        title(ax, [subj ' — missing']); continue
    end
    plot_tf_panel(ax, t_pull, fr, EPOCHS.(subj).tf_pull, coi_mask_pull, bwr, ...
        CAXIS_DB, FREQ_PLOT_MIN, FREQ_PLOT_MAX, tf_smooth, BANDS, subj, 'D', 0);
end
sgtitle('Per-Subject TF — Pull Onset (dB re rest)', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1600 800]);
print(fullfile(RESULTS_DIR,'TF_perSubject_Pull'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'TF_perSubject_Pull.fig'));

%% ===== SAVE =====
TF.subjects        = SUBJECTS;
TF.fr        = fr;
TF.fr         = fr;
TF.t_onset         = t_onset;
TF.t_pull          = t_pull;
TF.grand_tf_onset  = grand_tf_onset;
TF.grand_tf_pull   = grand_tf_pull;
TF.sem_tf_onset    = sem_tf_onset;
TF.sem_tf_pull     = sem_tf_pull;
TF.coi_mask_onset  = coi_mask_onset;
TF.coi_mask_pull   = coi_mask_pull;
TF.per_subj_onset  = all_tf_onset;    % [n_fr × n_samp × n_subj]
TF.per_subj_pull   = all_tf_pull;
TF.per_subj_baseline = all_bl;        % [n_fr × n_subj]
TF.params.freq_limits    = FREQ_LIMITS;
TF.params.vpo            = VOICES_PER_OCT;
TF.params.wavelet        = WAVELET_TYPE;
TF.params.baseline_win   = BASELINE_WIN;
TF.params.caxis_db       = CAXIS_DB;

save(fullfile(RESULTS_DIR, 'TF_allSubjects.mat'), 'TF', '-v7.3');
fprintf('\nSaved TF results to %s\n', RESULTS_DIR);

%% ===== LOCAL FUNCTION =====
% Must be placed after main script body (MATLAB R2016b+)

function plot_tf_panel(ax, t_vec, fr_vec, tf_mat, coi_mask, bwr_cmap, ...
    caxis_db, freq_min, freq_max, t_smooth, bands, title_str, ...
    event_label, event_t)

freq_mask = fr_vec >= freq_min & fr_vec <= freq_max;
fr_plot   = fr_vec(freq_mask);
tf_plot   = tf_mat(freq_mask, :);

% NaN-out COI region so pcolor renders it as background (white)
tf_masked          = tf_plot;
tf_masked(~coi_mask(freq_mask, :)) = NaN;

pcolor(ax, t_vec, fr_plot, t_smooth(tf_masked));
shading(ax, 'interp');
colormap(ax, bwr_cmap);
set(ax, 'YScale', 'log', 'FontSize', 10);
clim(ax, [-caxis_db, caxis_db]);
c = colorbar(ax); ylabel(c, 'dB re rest');
xline(ax, event_t, '--w', event_label, 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom');
xlabel(ax, 'Time (s)', 'FontSize', 11);
ylabel(ax, 'Frequency (Hz)', 'FontSize', 11);
title(ax, title_str, 'FontSize', 11);
ylim(ax, [freq_min, freq_max]);

for b = 1:numel(bands)
    if bands(b).freq >= freq_min && bands(b).freq <= freq_max
        yline(ax, bands(b).freq, ':', 'Color', [bands(b).col 0.6], ...
            'LineWidth', 1, 'Label', bands(b).name, ...
            'LabelVerticalAlignment', 'middle', 'FontSize', 8);
    end
end
end
