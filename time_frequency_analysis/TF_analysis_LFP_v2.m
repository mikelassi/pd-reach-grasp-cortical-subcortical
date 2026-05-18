%% ==========================================
%  LFP Time-Frequency Analysis  v2
%
%  Improvements over v1:
%   1. Baseline = rest_before ONLY (not rest_before + rest_after), so
%      post-movement beta rebound does not inflate the normalization
%      denominator.
%   2. dB normalization:  10*log10(power / baseline)  instead of
%      percent-change. dB is symmetric, additive, and better suited to
%      log-normally distributed power.
%   3. Single interpolation per trial: global median phase lengths are
%      computed across ALL trials and blocks in a first pass; each trial
%      is then interpolated ONCE directly to that common grid.  No further
%      temporal resampling is needed when averaging across blocks or
%      subjects.
%   4. pchip interpolation instead of linear, to avoid kinks at phase
%      boundaries and better preserve smooth power envelopes.
%   5. Trials are CONCATENATED into a subject-level pool before computing
%      statistics.  The median/mean is computed across all trials at once,
%      eliminating the mean-of-means bias that arises when blocks contain
%      different numbers of trials.
%   6. Fixed CWT frequency axis via cwtfilterbank with explicit
%      FrequencyLimits and VoicesPerOctave.  All blocks and subjects share
%      the same frequency vector, so no post-hoc truncation or
%      interpolation in the frequency dimension is needed.
%   7. All physical / analysis parameters declared at the top.
%   8. Colorbar limits set from a symmetric quantile (PCT_CAXIS), so
%      outlier trials do not compress the colour scale.
%   9. Non-normalised TF plot uses raw power units (not x100).
%
%  Author: Tommaso Marcantoni (revised by Michael Lassi)
%% ==========================================

close all; clear; clc;

%% ========== PARAMETERS  (edit here only) ==========

SUBJECTS        = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH       = 'I:\Parkinson_ReachGrasp\Reprocessing';      % root, no subject ID
PREPROC_SUB     = fullfile('Preprocessed', 'LFP');
RESULTS_DIR     = fullfile(BASE_PATH, 'RESULTS_final');

PHASE_NAMES     = {'Rest (pre)', 'Reach', 'Grasp', 'Pull', 'Rest (post)'};

% CWT  — analytic Morlet, extended to 80 Hz to capture gamma band
%   Reference paper used [1,50] Hz; here we extend to [1,80] Hz.
%   VoicesPerOctave = 12  →  ~76 log-spaced frequencies over [1,80] Hz
%     (≈ 6.3 octaves × 12 VPO; finer resolution than the paper's 9 VPO
%      so gamma bands are well-resolved without excessive computation).
%   COI contamination is NOT a concern above ~10 Hz: at 80 Hz the Morlet
%   is only ~125 ms wide, far shorter than any phase window.
%   NOTE: a 50 Hz notch filter was applied in preprocessing; expect a
%   small artefact notch at 50 Hz in the plots — this is expected.
FREQ_LIMITS     = [1, 80];    % Hz
VOICES_PER_OCT  = 12;         % ~76 log-spaced frequencies over [1,80] Hz
WAVELET_TYPE    = 'amor';     % analytic Morlet ('amor')

% Baseline
REST_OFFSET_S   = 0.5;       % seconds of pre-trial LFP included before phase 1
BASELINE_MODE   = 'pre';      % 'pre'  = rest-before only (PMBR visible, rest~0 dB)
                               % 'pre_post' = rest-before + rest-after (PMBR partially suppressed)

% Plotting
FREQ_PLOT_MIN   = 4;          % lower display limit (Hz)
                               % CWT is still computed from FREQ_LIMITS(1)=1 Hz for
                               % filter-bank stability, but below ~4 Hz the 10-cycle
                               % Morlet needs >2.5 s support — shorter than reach/grasp
                               % phases — so estimates are unreliable and noisy.
FREQ_PLOT_MAX   = 80;         % upper frequency limit for TF colour plots (Hz)
PCT_CAXIS       = 99;         % percentile for non-normalised TF colorbar clipping
CAXIS_DB_NORM   = 5;          % ± dB limit for normalised TF colormap
                               % Instantaneous CWT power has a wide range (individual beta
                               % bursts can hit ±10–15 dB), while phase-averaged modulation
                               % is ±1–3 dB.  A fixed ±3 dB limit shows physiological
                               % gradation; raise to 5 for subjects with very strong bursts.
TF_SMOOTH_SIGMA = [1, 8];     % 2-D Gaussian smoothing for TF display ONLY [freq, time] samples
                               % freq=1: minimal spectral blur; time=8: ~50-100 ms temporal
                               % smoothing to suppress episodic burst noise in single subjects.
                               % Set to [0, 0] to disable. Does NOT affect PSD or statistics.

%% ========== INIT ==========
n_subjects = numel(SUBJECTS);
n_phases   = numel(PHASE_NAMES);   % 5: Rest (pre), Reach, Grasp, Pull, Rest (post)

eeglab; close;   % initialise EEGLAB without GUI

%% ========== CROSS-SUBJECT STORAGE ==========
all_subj_wt             = cell(1, n_subjects);
all_subj_wt_no_norm     = cell(1, n_subjects);
all_subj_kin            = cell(1, n_subjects);
all_med_subj_baseline         = cell(1, n_subjects);
all_med_phases_psd_subj       = cell(n_phases, n_subjects);
all_med_phases_psd_norm_subj  = cell(n_phases, n_subjects);
all_med_psd_subj              = cell(1, n_subjects);
all_med_psd_norm_subj         = cell(1, n_subjects);
PSD_STRUCT = struct();

%% ========== PRE-PASS — universal phase grid (all subjects / trials) ==========
% global_med_len is computed ONCE from every trial in every block of every
% subject.  All trials are then warped to this single grid, so tot_length
% is identical across subjects and no second interpolation is needed for
% the grand average.
fprintf('Pre-pass: collecting phase lengths from all subjects...\n');
all_phase_lengths_global = [];
fs      = NaN;
fr_glob = [];
n_fr    = 0;

for s_pp = 1:n_subjects
    subj_pp    = SUBJECTS{s_pp};
    preproc_pp = fullfile(BASE_PATH, subj_pp, PREPROC_SUB);
    mat_pp     = dir(fullfile(preproc_pp, '*_LFP_trialsByRegionAndPhase.mat'));
    set_pp     = dir(fullfile(preproc_pp, '*.set'));

    if isnan(fs) && ~isempty(set_pp)
        tmp_lfp = pop_loadset('filename', set_pp(1).name, 'filepath', preproc_pp);
        fs      = tmp_lfp.srate;
        fb_tmp  = cwtfilterbank('SignalLength',  size(tmp_lfp.data, 2), ...
            'SamplingFrequency', fs, 'Wavelet', WAVELET_TYPE, ...
            'FrequencyLimits',   FREQ_LIMITS, 'VoicesPerOctave', VOICES_PER_OCT);
        fr_glob = centerFrequencies(fb_tmp);
        n_fr    = numel(fr_glob);
        clear tmp_lfp fb_tmp;
    end

    for f_pp = 1:numel(mat_pp)
        M_pp   = load(fullfile(preproc_pp, mat_pp(f_pp).name));
        LFP_ph = M_pp.LFP_phases;
        n_t_pp = size(LFP_ph, 1);
        blk    = zeros(n_t_pp, 5);
        for t_pp = 1:n_t_pp
            for p_pp = 1:5
                blk(t_pp, p_pp) = size(LFP_ph{t_pp, p_pp}, 2);
            end
        end
        all_phase_lengths_global = [all_phase_lengths_global; blk]; %#ok<AGROW>
    end
end

global_med_len   = round(median(all_phase_lengths_global, 1));  % 1×5, universal
tot_length       = sum(global_med_len);
rest_offset_samp = round(REST_OFFSET_S * fs);

% Percentage time axis: 0% = movement onset, 100% = movement offset.
% Rest-before appears at negative %, rest-after at >100%.
movement_samples = sum(global_med_len(2:4));
time_axis_pct    = ((1:tot_length) - 1 - global_med_len(1)) / movement_samples * 100;
t_onset_pct      = 0;
t_grasp_pct      = global_med_len(2)        / movement_samples * 100;
t_pull_pct       = sum(global_med_len(2:3)) / movement_samples * 100;
t_offset_pct     = 100;

fprintf('Universal grid (samples): %s  |  total = %d (%.2f s)\n', ...
    mat2str(global_med_len), tot_length, tot_length / fs);
fprintf('Phase markers: onset = 0%%  grasp = %.1f%%  pull = %.1f%%  offset = 100%%\n', ...
    t_grasp_pct, t_pull_pct);

%% ========== COLORMAPS (defined once) ==========
n_cmap = 256; half_c = n_cmap / 2;
ramp_c = linspace(0, 1, half_c)';
bwr_cmap = [[ramp_c, ramp_c, ones(half_c,1)]; ...
            [ones(half_c,1), flip(ramp_c), flip(ramp_c)]];  % blue–white–red

%% ========== SHARED PLOT CONSTANTS (defined once) ==========
col_rest      = [0.8 0.2 0.2];
col_reach     = [0.3 0.8 0.3];
col_grasp     = [0.2 0.6 1.0];
col_pull      = [0.6 0.1 0.9];
col_rest_post = [1.0 0.55 0.1];
phase_cols    = {col_rest, col_reach, col_grasp, col_pull, col_rest_post};

% Display-only 2-D Gaussian smoothing — sigma must be > 0
if all(TF_SMOOTH_SIGMA > 0)
    tf_smooth = @(M) imgaussfilt(double(M), TF_SMOOTH_SIGMA);
else
    tf_smooth = @(M) double(M);   % no-op when sigma = [0,0]
end

%% ========== SUBJECT LOOP ==========
for s = 1:n_subjects

    subject_id   = SUBJECTS{s};
    preproc_path = fullfile(BASE_PATH, subject_id, PREPROC_SUB);
    fprintf('\n===== Subject: %s =====\n', subject_id);

    mat_files = dir(fullfile(preproc_path, '*_LFP_trialsByRegionAndPhase.mat'));
    set_files = dir(fullfile(preproc_path, '*.set'));
    n_blocks  = numel(mat_files);




    % Contralateral hemisphere channel
    if ismember(subject_id, {'wue02', 'wue03'})
        ch = 2;
    else
        ch = 1;
    end

    % global_med_len, tot_length, rest_offset_samp, fr_glob, n_fr, fs
    % are all computed once in the Pre-Pass above.

    % ------------------------------------------------------------------
    % PASS — CWT, dB normalisation, single interpolation, concatenate
    % ------------------------------------------------------------------
    % Trial-pool arrays (3rd dimension grows as trials are appended)
    pool_wt          = zeros(n_fr, tot_length, 0);
    pool_wt_no_norm  = zeros(n_fr, tot_length, 0);
    pool_kin         = zeros(tot_length, 0);
    pool_psd         = zeros(n_fr, 0);
    pool_psd_norm    = zeros(n_fr, 0);
    pool_phase_psd      = zeros(n_phases, n_fr, 0);
    pool_phase_psd_norm = zeros(n_phases, n_fr, 0);
    pool_baseline    = zeros(n_fr, 0);

    for f = 1:n_blocks
        fprintf('  Block %d / %d\n', f, n_blocks);

        LFP = pop_loadset('filename', set_files(f).name, 'filepath', preproc_path);
        M   = load(fullfile(preproc_path, mat_files(f).name));

        baseline_LFP     = M.baseline_LFP;
        all_baseline_len = size(baseline_LFP, 2);
        LFP_trials       = M.LFP_trials;
        LFP_phases       = M.LFP_phases;
        KIN_trials       = M.kinematic_4LFP_trials;
        num_trials       = size(LFP_trials, 1);
        n_lfp_phases     = size(LFP_phases, 2);   % should be 5

        % CWT on the entire block — avoids edge effects at trial boundaries
        fb_block = cwtfilterbank('SignalLength',    size(LFP.data, 2), ...
            'SamplingFrequency', fs, ...
            'Wavelet',           WAVELET_TYPE, ...
            'FrequencyLimits',   FREQ_LIMITS, ...
            'VoicesPerOctave',   VOICES_PER_OCT);
        cwt_block     = wt(fb_block, LFP.data(ch, :));
        power_wt_glob = abs(cwt_block).^2;   % [n_fr x n_samples]

        % Sanity check: frequency axis must be consistent
        assert(numel(centerFrequencies(fb_block)) == n_fr, ...
            'Frequency axis length changed in block %d of subject %s', f, subject_id);

        start_idx = 0;   % tracks position in the block CWT (resets each block)

        for t = 1:num_trials
            kin_signal   = KIN_trials{t};

            % Buffers for this trial (interpolated to global grid)
            pow_wt_trial         = zeros(n_fr, tot_length);
            pow_wt_trial_no_norm = zeros(n_fr, tot_length);
            kin_trial_interp     = zeros(1, tot_length);

            phase_psd_trial      = zeros(n_phases, n_fr);
            phase_psd_norm_trial = zeros(n_phases, n_fr);
            mean_baseline_trial  = [];   % set when p == 1
            start_idx_trial      = NaN;  % start of phase 1 in block CWT
            raw_p1               = zeros(n_fr, 1);   % rest-before raw power
            raw_p5               = zeros(n_fr, 1);   % rest-after  raw power

            for p = 1:n_lfp_phases

                phase_length_nom = size(LFP_phases{t, p}, 2);

                % --- CWT start position ---
                if t == 1 && p == 1
                    start_idx = all_baseline_len - rest_offset_samp;
                else
                    start_idx = end_idx + 1;
                end

                % --- Kin start position ---
                if p == 1
                    start_kin = 1;
                else
                    start_kin = end_kin + 1;
                end

                % Clamp phase_length to whatever is available in BOTH signals.
                % A single clamp here prevents all downstream index errors.
                avail_pwr    = max(0, size(power_wt_glob, 2) - start_idx + 1);
                avail_kin    = max(0, length(kin_signal)     - start_kin  + 1);
                phase_length = min([phase_length_nom, avail_pwr, avail_kin]);

                end_idx = start_idx + phase_length - 1;
                end_kin = start_kin + phase_length - 1;

                % Target interpolated length for this phase
                tgt_len = global_med_len(p);

                if phase_length >= 2
                    pow_phase_no_norm = power_wt_glob(:, start_idx:end_idx);
                    if p == 1
                        % Geometric-mean baseline (mean in log domain).
                        % Guarantees rest PSD = 0 dB by construction:
                        %   mean_t[ log(P/geomean(P)) ] = 0 exactly.
                        % Arithmetic mean would give a systematic negative
                        % bias via Jensen's inequality because log is concave.
                        mean_baseline_trial = 10.^(mean(log10(max(pow_phase_no_norm, eps)), 2));
                        start_idx_trial     = start_idx;
                    end
                    pow_phase_db = 10 * log10(pow_phase_no_norm ./ mean_baseline_trial);

                    x_orig = (1:phase_length)';
                    x_new  = linspace(1, phase_length, tgt_len)';
                    pow_interp         = interp1(x_orig, pow_phase_db',      x_new, 'pchip', 'extrap')';
                    pow_interp_no_norm = interp1(x_orig, pow_phase_no_norm', x_new, 'pchip', 'extrap')';

                    kin_phase  = kin_signal(start_kin:end_kin);
                    kin_interp = interp1(x_orig, kin_phase(:), x_new, 'pchip', 'extrap');

                elseif phase_length == 1
                    pow_phase_no_norm = power_wt_glob(:, start_idx);
                    if p == 1
                        mean_baseline_trial = pow_phase_no_norm;
                        start_idx_trial     = start_idx;
                    end
                    pow_phase_db       = 10 * log10(pow_phase_no_norm ./ mean_baseline_trial);
                    pow_interp         = repmat(pow_phase_db,       1, tgt_len);
                    pow_interp_no_norm = repmat(pow_phase_no_norm,  1, tgt_len);
                    kin_interp         = repmat(kin_signal(start_kin), tgt_len, 1);

                else  % phase_length == 0: entirely out of bounds, fill with zeros
                    pow_phase_no_norm  = zeros(n_fr, 1);
                    pow_phase_db       = zeros(n_fr, 1);
                    if p == 1
                        mean_baseline_trial = ones(n_fr, 1);  % neutral fallback
                        start_idx_trial     = start_idx;
                    end
                    pow_interp         = zeros(n_fr, tgt_len);
                    pow_interp_no_norm = zeros(n_fr, tgt_len);
                    safe_idx = min(max(start_kin - 1, 1), length(kin_signal));
                    if isempty(kin_signal) || safe_idx < 1
                        last_kin = 0;
                    else
                        last_kin = kin_signal(safe_idx);
                    end
                    kin_interp = repmat(last_kin, tgt_len, 1);
                end

                % Cache raw power for modular baseline (used after the loop)
                if p == 1;             raw_p1 = pow_phase_no_norm; end
                if p == n_lfp_phases;  raw_p5 = pow_phase_no_norm; end

                if p <= n_phases
                    phase_psd_trial(p, :)      = mean(pow_phase_no_norm, 2);
                    phase_psd_norm_trial(p, :) = mean(pow_phase_db, 2);
                end

                % --- Assemble interpolated trial ---
                if p == 1
                    start_pi = 1;
                else
                    start_pi = end_pi + 1;
                end
                end_pi = start_pi + tgt_len - 1;

                pow_wt_trial(:, start_pi:end_pi)         = pow_interp;
                pow_wt_trial_no_norm(:, start_pi:end_pi) = pow_interp_no_norm;
                kin_trial_interp(start_pi:end_pi)        = kin_interp;

            end % phase loop

            % ------------------------------------------------------------------
            % Override baseline/normalisation for 'pre_post' mode.
            % raw_p1 and raw_p5 are the uninterpolated rest-before and
            % rest-after power spectra; their concatenation gives a baseline
            % anchored to a resting state that brackets the movement epoch.
            % Note: PMBR will appear less prominent compared to 'pre' mode
            % because high post-movement beta power inflates the denominator.
            % ------------------------------------------------------------------
            if strcmp(BASELINE_MODE, 'pre_post')
                mean_baseline_trial = 10.^(mean(log10(max(cat(2, raw_p1, raw_p5), eps)), 2));
                pow_wt_trial = 10 * log10(max(pow_wt_trial_no_norm, eps) ./ mean_baseline_trial);
                col_s = 1;
                for p_r = 1:n_phases
                    col_e = col_s + global_med_len(p_r) - 1;
                    phase_psd_norm_trial(p_r, :) = mean(pow_wt_trial(:, col_s:col_e), 2)';
                    col_s = col_e + 1;
                end
                psd_norm_trial = mean(pow_wt_trial, 2);
            end

            % Full-trial PSD: use the indices tracked by the phase loop so the
            % range is guaranteed to lie within power_wt_glob.  end_idx now
            % points to the last sample of phase 5.
            pow_trial_all  = power_wt_glob(:, start_idx_trial : end_idx);
            psd_trial      = mean(pow_trial_all, 2);
            if ~strcmp(BASELINE_MODE, 'pre_post')
                psd_norm_trial = mean(10 * log10(pow_trial_all ./ mean_baseline_trial), 2);
            end

            % Rescale kinematic signal to the displayed frequency range [FREQ_PLOT_MIN, FREQ_PLOT_MAX]
            kin_rng  = max(kin_trial_interp) - min(kin_trial_interp) + eps;
            kin_norm = (kin_trial_interp - min(kin_trial_interp)) / kin_rng;
            kin_norm = kin_norm * (FREQ_PLOT_MAX - FREQ_PLOT_MIN) + FREQ_PLOT_MIN;

            % --- Append trial to subject-level pool ---
            pool_wt(:, :, end+1)            = pow_wt_trial;
            pool_wt_no_norm(:, :, end+1)    = pow_wt_trial_no_norm;
            pool_kin(:, end+1)              = kin_norm(:);
            pool_psd(:, end+1)              = psd_trial;
            pool_psd_norm(:, end+1)         = psd_norm_trial;
            pool_phase_psd(:, :, end+1)     = phase_psd_trial;
            pool_phase_psd_norm(:, :, end+1)= phase_psd_norm_trial;
            pool_baseline(:, end+1)         = mean_baseline_trial;

            % Store per-trial results in PSD_STRUCT
            PSD_STRUCT.(subject_id).block(f).trial(t).psd_norm = psd_norm_trial;
            PSD_STRUCT.(subject_id).block(f).trial(t).psd      = psd_trial;
            for p = 1:n_phases
                PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd_norm = phase_psd_norm_trial(p, :)';
                PSD_STRUCT.(subject_id).block(f).trial(t).phase(p).psd      = phase_psd_trial(p, :)';
            end
        end % trial loop

        PSD_STRUCT.(subject_id).block(f).frequencies = fr_glob;

    end % block loop

    % ------------------------------------------------------------------
    % Subject-level statistics — computed across ALL pooled trials
    % (no mean-of-means; every trial contributes equally)
    % ------------------------------------------------------------------
    mean_wt_subj         = median(pool_wt, 3);       % median across trials — consistent with PSD aggregation
    mean_wt_subj_no_norm = median(pool_wt_no_norm, 3);
    mean_kin_subj        = median(pool_kin, 2);   % [tot_length x 1]

    for p = 1:n_phases
        all_med_phases_psd_subj{p, s}      = median(squeeze(pool_phase_psd(p, :, :)),      2);
        all_med_phases_psd_norm_subj{p, s} = median(squeeze(pool_phase_psd_norm(p, :, :)), 2);
    end
    all_med_psd_subj{s}      = median(pool_psd, 2);
    all_med_psd_norm_subj{s} = median(pool_psd_norm, 2);
    all_med_subj_baseline{s} = median(pool_baseline, 2);

    % Store for grand-average computation
    all_subj_wt{s}         = mean_wt_subj;
    all_subj_wt_no_norm{s} = mean_wt_subj_no_norm;
    all_subj_kin{s}        = mean_kin_subj;

    % ------------------------------------------------------------------
    % Per-subject TF plot (normalised)
    % ------------------------------------------------------------------
    freq_mask_s = fr_glob >= FREQ_PLOT_MIN & fr_glob <= FREQ_PLOT_MAX;
    wt_plot_s   = mean_wt_subj(freq_mask_s, :);
    freq_plot_s = fr_glob(freq_mask_s);

    figure('Name', ['Normalized TF - ' subject_id], 'NumberTitle', 'off');
    pcolor(time_axis_pct, freq_plot_s, tf_smooth(wt_plot_s));
    shading interp; colormap(bwr_cmap);
    set(gca, 'YScale', 'linear', 'FontSize', 12);
    xlabel('Normalized trial time (% of movement)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    caxis([-CAXIS_DB_NORM, CAXIS_DB_NORM]);
    c = colorbar; ylabel(c, 'Power (dB re rest)');
    xline(t_onset_pct,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_grasp_pct,  'Label', 'Grasp-start', 'Color', 'w', 'LineWidth', 2.5);
    xline(t_pull_pct,   'Label', 'Grasp-end',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_offset_pct, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
    title(['Normalized Power CWT - ' subject_id ' Ch' num2str(ch)], 'FontSize', 13);
    yyaxis right
    plot(time_axis_pct, mean_kin_subj, 'k', 'LineWidth', 2);
    ax_s = gca; ax_s.YAxis(2).Visible = 'off';

    out_dir_TF = fullfile(RESULTS_DIR, 'TF_pow_norm_phase_results_01');
    if ~exist(out_dir_TF, 'dir'), mkdir(out_dir_TF); end
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
    print(fullfile(out_dir_TF, sprintf('Averaged_TF_%s_all_blocks', subject_id)), '-dpng', '-r300');

    % ------------------------------------------------------------------
    % Per-subject TF plot — non-normalised (raw power)
    % ------------------------------------------------------------------
    wt_plot_s_nn = mean_wt_subj_no_norm(freq_mask_s, :);

    figure('Name', ['Non-Normalised TF - ' subject_id], 'NumberTitle', 'off');
    pcolor(time_axis_pct, freq_plot_s, tf_smooth(wt_plot_s_nn));
    shading interp; colormap(parula);
    set(gca, 'YScale', 'linear', 'FontSize', 12);
    xlabel('Normalized trial time (% of movement)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    clo_s_nn = quantile(wt_plot_s_nn(:), 1 - PCT_CAXIS/100);
    chi_s_nn = quantile(wt_plot_s_nn(:), PCT_CAXIS/100);
    caxis([clo_s_nn, chi_s_nn]);
    c = colorbar; ylabel(c, 'Power (a.u.)');
    xline(t_onset_pct,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_grasp_pct,  'Label', 'Grasp-start', 'Color', 'w', 'LineWidth', 2.5);
    xline(t_pull_pct,   'Label', 'Grasp-end',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_offset_pct, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
    title(['Non-Normalised Power CWT - ' subject_id ' Ch' num2str(ch)], 'FontSize', 13);
    yyaxis right
    plot(time_axis_pct, mean_kin_subj, 'k', 'LineWidth', 2);
    ax_s_nn = gca; ax_s_nn.YAxis(2).Visible = 'off';

    out_dir_TF_nn = fullfile(RESULTS_DIR, 'TF_pow_norm_phase_results_NO_NORM');
    if ~exist(out_dir_TF_nn, 'dir'), mkdir(out_dir_TF_nn); end
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
    print(fullfile(out_dir_TF_nn, sprintf('Averaged_TF_NoNorm_%s_all_blocks', subject_id)), '-dpng', '-r300');

    % ------------------------------------------------------------------
    % Per-subject PSD plot (normalised phases)
    % ------------------------------------------------------------------
    figure('Name', ['Phases Normalised PSD - ' subject_id], 'NumberTitle', 'off');
    hold on; grid on;
    psd_handles = gobjects(1, n_phases);
    for p = 1:n_phases
        psd_handles(p) = plot(fr_glob, all_med_phases_psd_norm_subj{p, s}, ...
            'Color', phase_cols{p}, 'LineWidth', 2.5, 'DisplayName', PHASE_NAMES{p});
    end
    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Power (dB re Rest)', 'FontSize', 12);
    xlim([FREQ_PLOT_MIN, FREQ_PLOT_MAX]);
    title(['Median Normalised PSD - Kinematic Phases - ' subject_id], 'FontSize', 14);
    legend(psd_handles, 'Location', 'best');
    set(gca, 'FontSize', 12);

    out_dir_norm = fullfile(RESULTS_DIR, 'Baseline_01_median', 'NORM');
    if ~exist(out_dir_norm, 'dir'), mkdir(out_dir_norm); end
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
    saveas(gcf, fullfile(out_dir_norm, ['PSD_norm_phases_' subject_id '.fig']));
    print(fullfile(out_dir_norm, ['PSD_norm_phases_' subject_id '.png']), '-dpng', '-r300');

end % subject loop

%% ========== SAVE PSD_STRUCT ==========
out_dir_struct = fullfile(RESULTS_DIR, 'Baseline_01_median');
if ~exist(out_dir_struct, 'dir'), mkdir(out_dir_struct); end
save(fullfile(out_dir_struct, 'PSD_STRUCT_allSubjects.mat'), 'PSD_STRUCT', '-v7.3');

%% ========== GRAND AVERAGE ACROSS SUBJECTS ==========
% All subjects share the universal grid from the Pre-Pass, so tot_length
% is identical for everyone — no second interpolation needed.
min_baseline_len = min(cellfun(@length, all_med_subj_baseline));

all_wt_norm    = cat(3, all_subj_wt{:});          % [n_fr × tot_length × n_subjects]
all_wt_no_norm = cat(3, all_subj_wt_no_norm{:});
all_kin_mat    = cat(2, all_subj_kin{:});          % [tot_length × n_subjects]

all_baseline_mat        = zeros(min_baseline_len, n_subjects);
all_phases_psd_mat      = zeros(n_phases, min_baseline_len, n_subjects);
all_phases_psd_norm_mat = zeros(n_phases, min_baseline_len, n_subjects);
for s = 1:n_subjects
    all_baseline_mat(:, s) = all_med_subj_baseline{s}(1:min_baseline_len);
    for p = 1:n_phases
        all_phases_psd_mat(p, :, s)      = all_med_phases_psd_subj{p, s}(1:min_baseline_len);
        all_phases_psd_norm_mat(p, :, s) = all_med_phases_psd_norm_subj{p, s}(1:min_baseline_len);
    end
end

grand_avg_wt    = mean(all_wt_norm,    3);
grand_avg_wt_nn = mean(all_wt_no_norm, 3);
grand_avg_kin   = mean(all_kin_mat,    2);
sem_wt          = std(all_wt_norm, 0, 3) / sqrt(n_subjects);

mean_phases_psd      = zeros(n_phases, min_baseline_len);
mean_phases_psd_norm = zeros(n_phases, min_baseline_len);
sem_phases_psd       = zeros(n_phases, min_baseline_len);
sem_phases_psd_norm  = zeros(n_phases, min_baseline_len);
for p = 1:n_phases
    vals      = squeeze(all_phases_psd_mat(p, :, :));       % [freq x subj]
    vals_norm = squeeze(all_phases_psd_norm_mat(p, :, :));
    mean_phases_psd(p, :)      = mean(vals,      2);
    mean_phases_psd_norm(p, :) = mean(vals_norm, 2);
    sem_phases_psd(p, :)       = std(vals,      0, 2) / sqrt(n_subjects);
    sem_phases_psd_norm(p, :)  = std(vals_norm, 0, 2) / sqrt(n_subjects);
end

% Phase boundaries are universal — reuse precomputed percentage values
% Frequency and kinematic overlays for grand-average plots
freq_mask_ga  = fr_glob >= FREQ_PLOT_MIN & fr_glob <= FREQ_PLOT_MAX;
freq_plot_ga  = fr_glob(freq_mask_ga);
kin_rng_ga    = max(grand_avg_kin) - min(grand_avg_kin) + eps;
kin_norm_ga   = (grand_avg_kin - min(grand_avg_kin)) / kin_rng_ga;
kin_norm_ga   = kin_norm_ga * (max(freq_plot_ga) - min(freq_plot_ga)) + min(freq_plot_ga);

out_dir_ga    = fullfile(RESULTS_DIR, 'TF_pow_norm_phase_results_01');
if ~exist(out_dir_ga, 'dir'), mkdir(out_dir_ga); end
out_dir_nn    = fullfile(RESULTS_DIR, 'TF_pow_norm_phase_results_NO_NORM');
if ~exist(out_dir_nn, 'dir'), mkdir(out_dir_nn); end
out_dir_psd   = fullfile(RESULTS_DIR, 'Baseline_01_median');
if ~exist(out_dir_psd, 'dir'), mkdir(out_dir_psd); end

%% ========== PLOT: Grand-Average Normalised TF ==========
figure('Name', 'Grand-Average LFP TF (dB, Rest-Normalised)', 'NumberTitle', 'off');
wt_plot_ga = grand_avg_wt(freq_mask_ga, :);
pcolor(time_axis_pct, freq_plot_ga, tf_smooth(wt_plot_ga));
shading interp; colormap(bwr_cmap);
set(gca, 'YScale', 'linear', 'FontSize', 12);
xlabel('Normalized trial time (% of movement)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
caxis([-CAXIS_DB_NORM, CAXIS_DB_NORM]);
c = colorbar; ylabel(c, 'Power (dB re Rest)');
xline(t_onset_pct,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
xline(t_grasp_pct,  'Label', 'Grasp Start', 'Color', 'w', 'LineWidth', 2.5);
xline(t_pull_pct,   'Label', 'Grasp End',   'Color', 'w', 'LineWidth', 2.5);
xline(t_offset_pct, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
title('Grand-Average LFP TF Power (dB, Rest-Normalised)', 'FontSize', 14);
yyaxis right
plot(time_axis_pct, kin_norm_ga, 'k', 'LineWidth', 2);
ax_ga = gca; ax_ga.YAxis(2).Visible = 'off';
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
print(fullfile(out_dir_ga, 'All_Subjects_Average_Power_TF_norm'), '-dpng', '-r300');
savefig(fullfile(out_dir_ga, 'All_Subjects_Average_Power_TF_norm.fig'));

%% ========== PLOT: Grand-Average Non-Normalised TF ==========
figure('Name', 'Grand-Average LFP TF (Not Normalised)', 'NumberTitle', 'off');
wt_plot_nn = grand_avg_wt_nn(freq_mask_ga, :);
pcolor(time_axis_pct, freq_plot_ga, tf_smooth(wt_plot_nn));
shading interp; colormap(parula);
set(gca, 'YScale', 'linear', 'FontSize', 12);
xlabel('Normalized trial time (% of movement)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
clo_nn = quantile(wt_plot_nn(:), 1 - PCT_CAXIS/100);
chi_nn = quantile(wt_plot_nn(:), PCT_CAXIS/100);
caxis([clo_nn, chi_nn]);
c = colorbar; ylabel(c, 'Power (a.u.)');
xline(t_onset_pct,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
xline(t_grasp_pct,  'Label', 'Grasp Start', 'Color', 'w', 'LineWidth', 2.5);
xline(t_pull_pct,   'Label', 'Grasp End',   'Color', 'w', 'LineWidth', 2.5);
xline(t_offset_pct, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
title('Grand-Average LFP TF Power (Not Normalised)', 'FontSize', 14);
yyaxis right
plot(time_axis_pct, kin_norm_ga, 'k', 'LineWidth', 2);
ax_nn = gca; ax_nn.YAxis(2).Visible = 'off';
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
print(fullfile(out_dir_nn, 'All_Subjects_Average_Power_TF_NO_NORM'), '-dpng', '-r300');
savefig(fullfile(out_dir_nn, 'All_Subjects_Average_Power_TF_NO_NORM.fig'));

%% ========== PLOT: Grand-Average PSD (Non-Normalised) ==========
common_fr_row = reshape(fr_glob(1:min_baseline_len), 1, []);   % ensure row vector

figure('Name', 'Grand-Average PSD (Non-Normalised)', 'Position', [100 100 1200 600]);
hold on; grid on;
for p = 1:n_phases
    mn = mean_phases_psd(p, :);
    se = sem_phases_psd(p, :);
    fill([common_fr_row, fliplr(common_fr_row)], [mn+se, fliplr(mn-se)], ...
        phase_cols{p}, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    plot(common_fr_row, mn, 'Color', phase_cols{p}, 'LineWidth', 2.5, 'DisplayName', PHASE_NAMES{p});
end
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (a.u.)', 'FontSize', 12);
xlim([FREQ_PLOT_MIN, FREQ_PLOT_MAX]);
title('Grand-Average LFP Marginal Power Spectral Density', 'FontSize', 14);
legend('Location', 'best');
set(gca, 'FontSize', 12);
print(fullfile(out_dir_psd, 'Mean_PSD_all_subjects'), '-dpng', '-r300');
saveas(gcf, fullfile(out_dir_psd, 'Mean_PSD_all_subjects.fig'));

%% ========== PLOT: Grand-Average PSD (Normalised, dB) ==========
figure('Name', 'Grand-Average PSD (dB, Rest-Normalised)', 'Position', [100 100 1200 600]);
hold on; grid on;
for p = 1:n_phases
    mn = mean_phases_psd_norm(p, :);
    se = sem_phases_psd_norm(p, :);
    fill([common_fr_row, fliplr(common_fr_row)], [mn+se, fliplr(mn-se)], ...
        phase_cols{p}, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    plot(common_fr_row, mn, 'Color', phase_cols{p}, 'LineWidth', 2.5, 'DisplayName', PHASE_NAMES{p});
end
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (dB re Rest)', 'FontSize', 12);
xlim([FREQ_PLOT_MIN, FREQ_PLOT_MAX]);
title('Grand-Average LFP Marginal Normalised PSD (dB re Rest)', 'FontSize', 14);
legend('Location', 'best');
set(gca, 'FontSize', 12);
print(fullfile(out_dir_psd, 'Mean_PSD_norm_all_subjects'), '-dpng', '-r300');
saveas(gcf, fullfile(out_dir_psd, 'Mean_PSD_norm_all_subjects.fig'));