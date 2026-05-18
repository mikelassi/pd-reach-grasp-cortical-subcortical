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
%  Author: Tommaso Marcantoni (revised)
%% ==========================================

close all; clear; clc;

%% ========== PARAMETERS  (edit here only) ==========

SUBJECTS        = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH       = 'H:\Parkinson_ReachGrasp\Reprocessing';      % root, no subject ID
PREPROC_SUB     = fullfile('Preprocessed', 'LFP');
RESULTS_DIR     = fullfile(BASE_PATH, 'RESULTS_final');

PHASE_NAMES     = {'Rest', 'Reach', 'Grasp', 'Pull'};

% CWT
FREQ_LIMITS     = [1, 200];   % Hz  (lower limit must be reachable for shortest block)
VOICES_PER_OCT  = 12;         % controls frequency resolution

% Baseline
REST_OFFSET_S   = 0.25;       % seconds of pre-trial LFP included before phase 1

% Plotting
FREQ_PLOT_MAX   = 90;         % upper frequency limit for TF colour plots (Hz)
PCT_CAXIS       = 99;         % symmetric percentile for colorbar clipping

%% ========== INIT ==========
n_subjects = numel(SUBJECTS);
n_phases   = numel(PHASE_NAMES);   % 4: Rest, Reach, Grasp, Pull

eeglab; close;   % initialise EEGLAB without GUI

%% ========== CROSS-SUBJECT STORAGE ==========
all_subj_wt             = cell(1, n_subjects);
all_subj_wt_no_norm     = cell(1, n_subjects);
all_subj_kin            = cell(1, n_subjects);
all_subj_fr             = cell(1, n_subjects);
all_kin_latencies       = zeros(4, n_subjects);
all_med_subj_baseline         = cell(1, n_subjects);
all_med_phases_psd_subj       = cell(n_phases, n_subjects);
all_med_phases_psd_norm_subj  = cell(n_phases, n_subjects);
all_med_psd_subj              = cell(1, n_subjects);
all_med_psd_norm_subj         = cell(1, n_subjects);
PSD_STRUCT = struct();

grand_fs = NaN;   % filled from first subject; assumed equal across subjects

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

    % ------------------------------------------------------------------
    % PASS 1 — collect phase lengths and sampling rate
    % ------------------------------------------------------------------
    all_phase_lengths = [];   % (total_trials_across_blocks) x 5
    fs = [];

    for f = 1:n_blocks
        if isempty(fs)
            LFP_tmp = pop_loadset('filename', set_files(f).name, 'filepath', preproc_path);
            fs = LFP_tmp.srate;
            clear LFP_tmp;
        end
        M_tmp = load(fullfile(preproc_path, mat_files(f).name));
        LFP_ph = M_tmp.LFP_phases;
        n_t = size(LFP_ph, 1);
        block_len = zeros(n_t, 5);
        for t = 1:n_t
            for p = 1:5
                block_len(t, p) = size(LFP_ph{t, p}, 2);
            end
        end
        all_phase_lengths = [all_phase_lengths; block_len]; %#ok<AGROW>
    end

    % Global median phase lengths (computed across ALL trials in this subject)
    global_med_len = round(median(all_phase_lengths, 1));   % 1 x 5
    tot_length     = sum(global_med_len);
    rest_offset_samp = round(REST_OFFSET_S * fs);

    % Phase boundary times (seconds) — used for plot markers
    t_onset  = global_med_len(1) / fs;
    t_grasp  = t_onset  + global_med_len(2) / fs;
    t_pull   = t_grasp  + global_med_len(3) / fs;
    t_offset = t_pull   + global_med_len(4) / fs;

    % ------------------------------------------------------------------
    % Determine fixed CWT frequency axis via cwtfilterbank.
    % With explicit FrequencyLimits and VoicesPerOctave the center
    % frequencies depend only on those two parameters, not on SignalLength,
    % so all blocks share the same frequency vector.
    % ------------------------------------------------------------------
    LFP_ref = pop_loadset('filename', set_files(1).name, 'filepath', preproc_path);
    fb_ref  = cwtfilterbank('SignalLength',    size(LFP_ref.data, 2), ...
                            'SamplingFrequency', fs, ...
                            'Wavelet',           'morse', ...
                            'FrequencyLimits',   FREQ_LIMITS, ...
                            'VoicesPerOctave',   VOICES_PER_OCT);
    fr_glob = centerFrequencies(fb_ref);
    n_fr    = numel(fr_glob);
    clear LFP_ref fb_ref;

    if isnan(grand_fs)
        grand_fs = fs;
    end

    % ------------------------------------------------------------------
    % PASS 2 — CWT, dB normalisation, single interpolation, concatenate
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
                                 'Wavelet',           'morse', ...
                                 'FrequencyLimits',   FREQ_LIMITS, ...
                                 'VoicesPerOctave',   VOICES_PER_OCT);
        cwt_block     = wt(fb_block, LFP.data(ch, :));
        power_wt_glob = abs(cwt_block).^2;   % [n_fr x n_samples]

        % Sanity check: frequency axis must be consistent
        assert(numel(centerFrequencies(fb_block)) == n_fr, ...
            'Frequency axis length changed in block %d of subject %s', f, subject_id);

        start_idx = 0;   % tracks position in the block CWT (resets each block)

        for t = 1:num_trials
            trial_data   = LFP_trials{t}(ch, :);
            kin_signal   = KIN_trials{t};
            trial_length = size(trial_data, 2);

            % Buffers for this trial (interpolated to global grid)
            pow_wt_trial         = zeros(n_fr, tot_length);
            pow_wt_trial_no_norm = zeros(n_fr, tot_length);
            kin_trial_interp     = zeros(1, tot_length);

            phase_psd_trial      = zeros(n_phases, n_fr);
            phase_psd_norm_trial = zeros(n_phases, n_fr);
            mean_baseline_trial  = [];   % set when p == 1

            for p = 1:n_lfp_phases

                phase_length = size(LFP_phases{t, p}, 2);

                % Index range of this phase in the block CWT
                if t == 1 && p == 1
                    % First trial: overlap rest_offset_samp with the
                    % pre-trial baseline period
                    start_idx = all_baseline_len - rest_offset_samp;
                else
                    start_idx = end_idx + 1;
                end
                end_idx = start_idx + phase_length - 1;

                pow_phase_no_norm = power_wt_glob(:, start_idx:end_idx);

                % --- Baseline: rest_before ONLY (p == 1) ---
                % Using only pre-movement rest ensures post-movement beta
                % rebound does not inflate the normalisation denominator.
                if p == 1
                    mean_baseline_trial = mean(pow_phase_no_norm, 2);   % [n_fr x 1]
                end

                % dB normalisation against rest_before
                % 10*log10 is symmetric around 0, unbounded in both directions,
                % and better suited to log-normally distributed power than
                % percent-change.
                pow_phase_db = 10 * log10(pow_phase_no_norm ./ mean_baseline_trial);

                % Store phase PSD (phases 1-4; phase 5 = rest_after is not
                % stored separately since we use rest_before for baseline)
                if p <= n_phases
                    phase_psd_trial(p, :)      = mean(pow_phase_no_norm, 2);
                    phase_psd_norm_trial(p, :) = mean(pow_phase_db, 2);
                    % Note: phase 1 (rest_before) normalised against its own
                    % mean gives ~0 dB by construction — confirms flat baseline.
                end

                % Full-trial PSD (stored only once, when p == 1, covering the
                % whole trial from start_idx to start_idx+trial_length-1)
                if p == 1
                    pow_trial_all  = power_wt_glob(:, start_idx : start_idx + trial_length - 1);
                    psd_trial      = mean(pow_trial_all, 2);
                    psd_norm_trial = mean(10 * log10(pow_trial_all ./ mean_baseline_trial), 2);
                end

                % --- Single interpolation to global median phase length ---
                tgt_len = global_med_len(p);
                x_orig  = (1 : phase_length)';
                x_new   = linspace(1, phase_length, tgt_len)';

                % pchip preserves smooth power envelopes without kinks at
                % phase boundaries and without the overshoot of spline.
                pow_interp         = interp1(x_orig, pow_phase_db',         x_new, 'pchip', 'extrap')';
                pow_interp_no_norm = interp1(x_orig, pow_phase_no_norm',    x_new, 'pchip', 'extrap')';

                % --- Kinematic signal for this phase ---
                if p == 1
                    start_kin = 1;
                else
                    start_kin = end_kin + 1;
                end
                end_kin   = start_kin + phase_length - 1;
                kin_phase = kin_signal(start_kin : end_kin);
                kin_interp = interp1(x_orig, kin_phase(:), x_new, 'pchip', 'extrap');

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

            % Rescale kinematic signal to the displayed frequency range
            kin_rng  = max(kin_trial_interp) - min(kin_trial_interp) + eps;
            kin_norm = (kin_trial_interp - min(kin_trial_interp)) / kin_rng;
            kin_norm = kin_norm * (max(fr_glob) - min(fr_glob)) + min(fr_glob);

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
    mean_wt_subj         = mean(pool_wt, 3);
    mean_wt_subj_no_norm = mean(pool_wt_no_norm, 3);
    mean_kin_subj        = mean(pool_kin, 2);   % [tot_length x 1]

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
    all_subj_fr{s}         = fr_glob;

    % Kinematic latencies (from global median lengths, same for all blocks)
    all_kin_latencies(1, s) = t_onset;
    all_kin_latencies(2, s) = t_grasp;
    all_kin_latencies(3, s) = t_pull;
    all_kin_latencies(4, s) = t_offset;

    % ------------------------------------------------------------------
    % Per-subject TF plot
    % ------------------------------------------------------------------
    time_axis_subj = (1:tot_length) / fs;
    freq_mask_s    = fr_glob <= FREQ_PLOT_MAX;
    wt_plot_s      = mean_wt_subj(freq_mask_s, :);
    freq_plot_s    = fr_glob(freq_mask_s);

    figure('Name', ['Normalized TF - ' subject_id], 'NumberTitle', 'off');
    pcolor(time_axis_subj, freq_plot_s, wt_plot_s);
    shading interp; colormap jet;
    set(gca, 'YScale', 'linear', 'FontSize', 12);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    clo_s = quantile(wt_plot_s(:), 1 - PCT_CAXIS/100);
    chi_s = quantile(wt_plot_s(:), PCT_CAXIS/100);
    caxis([clo_s, chi_s]);
    c = colorbar; ylabel(c, 'Power (dB re rest)');
    xline(t_onset,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_grasp,  'Label', 'Grasp-start', 'Color', 'w', 'LineWidth', 2.5);
    xline(t_pull,   'Label', 'Grasp-end',   'Color', 'w', 'LineWidth', 2.5);
    xline(t_offset, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
    title(['Normalized Power CWT - ' subject_id ' Ch' num2str(ch)], 'FontSize', 13);
    yyaxis right
    plot(time_axis_subj, mean_kin_subj, 'k', 'LineWidth', 2);
    ax_s = gca; ax_s.YAxis(2).Visible = 'off';

    out_dir_TF = fullfile(RESULTS_DIR, 'TF_pow_norm_phase_results_01');
    if ~exist(out_dir_TF, 'dir'), mkdir(out_dir_TF); end
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
    print(fullfile(out_dir_TF, sprintf('Averaged_TF_%s_all_blocks', subject_id)), '-dpng', '-r300');

    % ------------------------------------------------------------------
    % Per-subject PSD plot (normalised phases)
    % ------------------------------------------------------------------
    col_rest  = [0.8 0.2 0.2];
    col_reach = [0.3 0.8 0.3];
    col_grasp = [0.2 0.6 1.0];
    col_pull  = [0.6 0.1 0.9];
    phase_cols = {col_rest, col_reach, col_grasp, col_pull};

    figure('Name', ['Phases Normalised PSD - ' subject_id], 'NumberTitle', 'off');
    hold on; grid on;
    psd_handles = gobjects(1, n_phases);
    for p = 1:n_phases
        psd_handles(p) = plot(fr_glob, all_med_phases_psd_norm_subj{p, s}, ...
            'Color', phase_cols{p}, 'LineWidth', 2.5, 'DisplayName', PHASE_NAMES{p});
    end
    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Power (dB re Rest)', 'FontSize', 12);
    xlim([min(fr_glob), 50]);
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
% Since all subjects now share the same frequency axis (fixed via
% cwtfilterbank parameters), min_freqs should equal n_fr; the truncation
% is kept as a safeguard.
min_freqs        = min(cellfun(@(x) size(x, 1), all_subj_wt));
median_len_ga    = round(median(cellfun(@(x) size(x, 2), all_subj_wt)));
min_baseline_len = min(cellfun(@length, all_med_subj_baseline));

common_fr    = all_subj_fr{1}(1:min_freqs);   % reference frequency axis
time_axis_ga = (1 : median_len_ga) / grand_fs;

% 3-D matrices: [freq x time x subjects] and [freq x subjects]
all_wt_norm   = zeros(min_freqs, median_len_ga, n_subjects);
all_wt_no_norm= zeros(min_freqs, median_len_ga, n_subjects);
all_kin_mat   = zeros(median_len_ga, n_subjects);
all_baseline_mat         = zeros(min_baseline_len, n_subjects);
all_phases_psd_mat       = zeros(n_phases, min_baseline_len, n_subjects);
all_phases_psd_norm_mat  = zeros(n_phases, min_baseline_len, n_subjects);

for s = 1:n_subjects
    wt_s    = all_subj_wt{s}(1:min_freqs, :);
    wt_s_nn = all_subj_wt_no_norm{s}(1:min_freqs, :);
    kin_s   = all_subj_kin{s};

    or_len   = size(wt_s, 2);
    t_orig   = 1:or_len;
    t_target = linspace(1, or_len, median_len_ga);

    all_wt_norm(:, :, s)   = interp1(t_orig, wt_s',    t_target, 'pchip', 'extrap')';
    all_wt_no_norm(:, :, s)= interp1(t_orig, wt_s_nn', t_target, 'pchip', 'extrap')';
    all_kin_mat(:, s)      = interp1(t_orig, kin_s(:)', t_target, 'pchip', 'extrap');

    all_baseline_mat(:, s) = all_med_subj_baseline{s}(1:min_baseline_len);
    for p = 1:n_phases
        all_phases_psd_mat(p, :, s)      = all_med_phases_psd_subj{p, s}(1:min_baseline_len);
        all_phases_psd_norm_mat(p, :, s) = all_med_phases_psd_norm_subj{p, s}(1:min_baseline_len);
    end
end

grand_avg_wt   = mean(all_wt_norm,    3);
grand_avg_wt_nn= mean(all_wt_no_norm, 3);
grand_avg_kin  = mean(all_kin_mat,    2);
sem_wt         = std(all_wt_norm, 0, 3) / sqrt(n_subjects);

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

% Mean kinematic event latencies across subjects
mean_onset_all  = mean(all_kin_latencies(1, :));
mean_grasp_all  = mean(all_kin_latencies(2, :));
mean_pull_all   = mean(all_kin_latencies(3, :));
mean_offset_all = mean(all_kin_latencies(4, :));

% Frequency and kinematic overlays for grand-average plots
freq_mask_ga  = common_fr <= FREQ_PLOT_MAX;
freq_plot_ga  = common_fr(freq_mask_ga);
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
pcolor(time_axis_ga, freq_plot_ga, wt_plot_ga);
shading interp; colormap jet;
set(gca, 'YScale', 'linear', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
clo_ga = quantile(wt_plot_ga(:), 1 - PCT_CAXIS/100);
chi_ga = quantile(wt_plot_ga(:), PCT_CAXIS/100);
caxis([clo_ga, chi_ga]);
c = colorbar; ylabel(c, 'Power (dB re Rest)');
xline(mean_onset_all,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
xline(mean_grasp_all,  'Label', 'Grasp Start', 'Color', 'w', 'LineWidth', 2.5);
xline(mean_pull_all,   'Label', 'Grasp End',   'Color', 'w', 'LineWidth', 2.5);
xline(mean_offset_all, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
title('Grand-Average LFP TF Power (dB, Rest-Normalised)', 'FontSize', 14);
yyaxis right
plot(time_axis_ga, kin_norm_ga, 'k', 'LineWidth', 2);
ax_ga = gca; ax_ga.YAxis(2).Visible = 'off';
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
print(fullfile(out_dir_ga, 'All_Subjects_Average_Power_TF_norm'), '-dpng', '-r300');
savefig(fullfile(out_dir_ga, 'All_Subjects_Average_Power_TF_norm.fig'));

%% ========== PLOT: Grand-Average Non-Normalised TF ==========
figure('Name', 'Grand-Average LFP TF (Not Normalised)', 'NumberTitle', 'off');
wt_plot_nn = grand_avg_wt_nn(freq_mask_ga, :);
pcolor(time_axis_ga, freq_plot_ga, wt_plot_nn);
shading interp; colormap jet;
set(gca, 'YScale', 'linear', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
clo_nn = quantile(wt_plot_nn(:), 1 - PCT_CAXIS/100);
chi_nn = quantile(wt_plot_nn(:), PCT_CAXIS/100);
caxis([clo_nn, chi_nn]);
c = colorbar; ylabel(c, 'Power (a.u.)');
xline(mean_onset_all,  'Label', 'Mov-ONSET',   'Color', 'w', 'LineWidth', 2.5);
xline(mean_grasp_all,  'Label', 'Grasp Start', 'Color', 'w', 'LineWidth', 2.5);
xline(mean_pull_all,   'Label', 'Grasp End',   'Color', 'w', 'LineWidth', 2.5);
xline(mean_offset_all, 'Label', 'Mov-OFFSET',  'Color', 'w', 'LineWidth', 2.5);
title('Grand-Average LFP TF Power (Not Normalised)', 'FontSize', 14);
yyaxis right
plot(time_axis_ga, kin_norm_ga, 'k', 'LineWidth', 2);
ax_nn = gca; ax_nn.YAxis(2).Visible = 'off';
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
print(fullfile(out_dir_nn, 'All_Subjects_Average_Power_TF_NO_NORM'), '-dpng', '-r300');
savefig(fullfile(out_dir_nn, 'All_Subjects_Average_Power_TF_NO_NORM.fig'));

%% ========== PLOT: Grand-Average PSD (Non-Normalised) ==========
col_rest  = [0.8 0.2 0.2];
col_reach = [0.3 0.8 0.3];
col_grasp = [0.2 0.6 1.0];
col_pull  = [0.6 0.1 0.9];
phase_cols = {col_rest, col_reach, col_grasp, col_pull};
common_fr_row = common_fr(:)';

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
xlim([0, 80]);
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
xlim([0, 80]);
title('Grand-Average LFP Marginal Normalised PSD (dB re Rest)', 'FontSize', 14);
legend('Location', 'best');
set(gca, 'FontSize', 12);
print(fullfile(out_dir_psd, 'Mean_PSD_norm_all_subjects'), '-dpng', '-r300');
saveas(gcf, fullfile(out_dir_psd, 'Mean_PSD_norm_all_subjects.fig'));
