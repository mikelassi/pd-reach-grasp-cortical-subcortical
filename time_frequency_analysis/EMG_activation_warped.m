%% EMG_activation_warped.m
% =========================================================================
% Time-warped average EMG activation across the full movement cycle.
%
% Mirrors TF_analysis_LFP_v2.m exactly in how trials are handled:
%   1. Pre-pass: collect all phase lengths from *_LFP_trialsByRegionAndPhase.mat
%      across all subjects/blocks → global median phase grid.
%   2. Subject loop: for each block, preprocess the full-block EMG envelope,
%      then extract phase segments using the SAME sequential index tracking
%      used in v2. Each phase is pchip-resampled to the global median length.
%   3. Trials are pooled (not block-averaged) before computing the median
%      activation per subject, eliminating mean-of-means bias.
%   4. Grand average ± SEM across subjects.
%
% Phases (from *_LFP_trialsByRegionAndPhase.mat):
%   1 — Rest (pre)   : window before movement onset
%   2 — Reach        : A → C
%   3 — Grasp        : C → D
%   4 — Pull         : D → F
%   5 — Rest (post)  : window after movement offset
%
% Outputs (RESULTS_DIR/EMG_Activation_Warped/):
%   EMG_grandAvg_warped.png / .fig
%   EMG_perSubject_warped.png / .fig
%   EMG_activation_warped.mat
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
SUBJECTS    = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
LFP_SUB     = fullfile('Preprocessed', 'LFP');      % *_LFP_trialsByRegionAndPhase.mat live here
EMG_SUB     = fullfile('01_Extracted', 'EMG_KIN');  % raw EMG signal
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'EMG_Activation_Warped');

PHASE_NAMES = {'Rest (pre)', 'Reach', 'Grasp', 'Pull', 'Rest (post)'};
n_phases    = numel(PHASE_NAMES);

% REST_OFFSET_S: seconds of pre-trial signal included as phase 1.
% Must match the value used in TF_analysis_LFP_v2.m so that the grid
% produced here has the same phase-1 length and time axis.
REST_OFFSET_S = 0.5;

% Target muscles — partial label matches, case-insensitive
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};
TARGET_MUSCLES    = {'iod', 'triceps', 'deltoid'};
MUSCLE_LABELS     = {'IOD', 'Triceps', 'Deltoid'};
n_muscles         = numel(TARGET_MUSCLES);

% EMG preprocessing
BP_LOW  = 20;  BP_HIGH = 180;  BP_ORDER = 4;
NOTCH_FREQS = [50, 100];  NOTCH_BW = 2;  NOTCH_ORDER = 2;
ENV_CUTOFF  = 5;  ENV_ORDER = 4;
NORM_PCTILE = 95;

% Display
MUSCLE_COLS = [0.15 0.65 0.25;   % IOD     — green
               0.80 0.30 0.10;   % Triceps — red-orange
               0.20 0.40 0.80];  % Deltoid — blue
COL_SHADE   = 0.25;
PHASE_COLS  = {[0.8 0.2 0.2],[0.3 0.8 0.3],[0.2 0.6 1.0],[0.6 0.1 0.9],[1.0 0.55 0.1]};

%% ===== BIOSIG STUB FIX =====
if ~isempty(which('eeglab'))
    bs = fullfile(fileparts(which('eeglab')), ...
        'plugins', 'Biosig3.8.4', 'biosig', 'maybe-missing');
    if exist(bs, 'dir'), rmpath(bs); end
end

%% ===== INIT =====
eeglab nogui;
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
n_subj = numel(SUBJECTS);

%% ===== PRE-PASS: global median phase lengths =====
% Identical to TF_analysis_LFP_v2.m pre-pass — reads the .mat files only
% for their phase lengths; no signal processing.
fprintf('Pre-pass: collecting phase lengths from all subjects...\n');
all_phase_lengths_global = [];
fs_lfp = NaN;

for s_pp = 1:n_subj
    lfp_path_pp = fullfile(BASE_PATH, SUBJECTS{s_pp}, LFP_SUB);
    mat_pp = dir(fullfile(lfp_path_pp, '*_LFP_trialsByRegionAndPhase.mat'));
    set_pp = dir(fullfile(lfp_path_pp, '*.set'));

    if isnan(fs_lfp) && ~isempty(set_pp)
        tmp = pop_loadset('filename', set_pp(1).name, 'filepath', lfp_path_pp);
        fs_lfp = tmp.srate;
        clear tmp;
    end

    for f_pp = 1:numel(mat_pp)
        M_pp   = load(fullfile(lfp_path_pp, mat_pp(f_pp).name));
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

global_med_len   = round(median(all_phase_lengths_global, 1));  % [1×5]
tot_length       = sum(global_med_len);
rest_offset_samp_lfp = round(REST_OFFSET_S * fs_lfp);

% Percentage time axis (0% = movement onset, 100% = movement offset)
movement_samples = sum(global_med_len(2:4));
time_axis_pct    = ((1:tot_length) - 1 - global_med_len(1)) / movement_samples * 100;
t_grasp_pct      = global_med_len(2)        / movement_samples * 100;
t_pull_pct       = sum(global_med_len(2:3)) / movement_samples * 100;

fprintf('Global grid (LFP samp): %s  →  total = %d samp  (%.2f s at %d Hz)\n', ...
    mat2str(global_med_len), tot_length, tot_length / fs_lfp, fs_lfp);

%% ===== PER-SUBJECT STORAGE =====
subj_median_emg = cell(n_muscles, n_subj);   % each: [1 × tot_length]

%% ===== SUBJECT LOOP =====
for s = 1:n_subj
    subj     = SUBJECTS{s};
    lfp_path = fullfile(BASE_PATH, subj, LFP_SUB);
    emg_path = fullfile(BASE_PATH, subj, EMG_SUB);

    mat_files = dir(fullfile(lfp_path, '*_LFP_trialsByRegionAndPhase.mat'));
    emg_files = dir(fullfile(emg_path, '*.set'));

    if isempty(mat_files)
        warning('Subject %s: no LFP .mat files — skipping.', subj);
        continue;
    end
    if isempty(emg_files)
        warning('Subject %s: no EMG .set files — skipping.', subj);
        continue;
    end

    fprintf('\n===== Subject: %s  (%d block(s)) =====\n', subj, numel(mat_files));
    use_dominant = ismember(subj, DOMINANT_SUBJECTS);

    % Trial pool: [n_muscles × tot_length × n_trials_total]
    pool_emg = zeros(n_muscles, tot_length, 0);

    for f = 1:min(numel(mat_files), numel(emg_files))
        fprintf('  Block %d\n', f);

        % ---- Load LFP .mat for phase timing ----
        M         = load(fullfile(lfp_path, mat_files(f).name));
        LFP_phases = M.LFP_phases;
        num_trials = size(LFP_phases, 1);
        all_baseline_len = size(M.baseline_LFP, 2);

        % ---- Load and preprocess full-block EMG ----
        EMG    = pop_loadset('filename', emg_files(f).name, 'filepath', emg_path);
        fs_emg = EMG.srate;

        % Channel selection
        ch_labels = lower({EMG.chanlocs.labels});
        if use_dominant
            sel_mask = contains(ch_labels, '_dominant_emg');
        else
            sel_mask = contains(ch_labels, '_emg');
        end
        ch_all = find(sel_mask);

        muscle_ch = zeros(1, n_muscles);
        for m = 1:n_muscles
            hits = find(contains(ch_labels(ch_all), TARGET_MUSCLES{m}));
            if ~isempty(hits), muscle_ch(m) = ch_all(hits(1)); end
        end
        if ~any(muscle_ch)
            warning('  Block %d: no target EMG channels — skipping.', f);
            continue;
        end

        % Build filters at EMG fs
        [b_bp,  a_bp]  = butter(BP_ORDER,  [BP_LOW, BP_HIGH]/(fs_emg/2), 'bandpass');
        [b_env, a_env]  = butter(ENV_ORDER, ENV_CUTOFF/(fs_emg/2),        'low');
        notch_b = cell(numel(NOTCH_FREQS), 2);
        for ni = 1:numel(NOTCH_FREQS)
            fn = NOTCH_FREQS(ni);
            if fn < fs_emg/2
                bw = NOTCH_BW/(fs_emg/2); fc = fn/(fs_emg/2);
                [notch_b{ni,1}, notch_b{ni,2}] = butter(NOTCH_ORDER, [fc-bw/2, fc+bw/2], 'stop');
            end
        end

        % Process full-block envelope for each muscle [n_muscles × n_samp_emg]
        env_block = zeros(n_muscles, EMG.pnts);
        for m = 1:n_muscles
            if muscle_ch(m) == 0, continue; end
            sig = double(EMG.data(muscle_ch(m),:));
            sig = sig - mean(sig);
            sig = filtfilt(b_bp, a_bp, sig);
            for ni = 1:numel(NOTCH_FREQS)
                if ~isempty(notch_b{ni,1})
                    sig = filtfilt(notch_b{ni,1}, notch_b{ni,2}, sig);
                end
            end
            sig = abs(sig);
            sig = filtfilt(b_env, a_env, sig);
            p95 = prctile(sig, NORM_PCTILE);
            if p95 > 0, sig = sig / p95; end
            env_block(m, :) = sig;
        end

        % fs ratio: scale LFP sample indices to EMG sample space
        fs_ratio = fs_emg / fs_lfp;

        % ---- TRIAL LOOP (mirrors TF_analysis_LFP_v2.m exactly) ----
        start_idx_lfp = 0;   % LFP sample counter (reset each block)

        for t = 1:num_trials
            emg_trial = zeros(n_muscles, tot_length);

            for p = 1:n_phases
                phase_len_lfp_nom = size(LFP_phases{t, p}, 2);

                % Advance LFP index (same logic as TF v2)
                if t == 1 && p == 1
                    start_idx_lfp = all_baseline_len - rest_offset_samp_lfp;
                else
                    start_idx_lfp = end_idx_lfp + 1;
                end

                % Convert LFP boundaries to EMG sample space
                start_emg = round(start_idx_lfp  * fs_ratio) + 1;   % 1-based
                end_emg_nom = round((start_idx_lfp + phase_len_lfp_nom - 1) * fs_ratio) + 1;

                % Clamp to valid EMG range
                start_emg   = max(1, min(start_emg,   EMG.pnts));
                end_emg_nom = max(1, min(end_emg_nom, EMG.pnts));
                phase_len_emg = end_emg_nom - start_emg + 1;

                tgt_len = global_med_len(p);

                if phase_len_emg >= 2
                    seg = env_block(:, start_emg : start_emg + phase_len_emg - 1);  % [n_muscles × len]
                    x_orig = (1:phase_len_emg)';
                    x_new  = linspace(1, phase_len_emg, tgt_len)';
                    seg_warp = interp1(x_orig, seg', x_new, 'pchip', 'extrap')';  % [n_muscles × tgt]
                elseif phase_len_emg == 1
                    seg_warp = repmat(env_block(:, start_emg), 1, tgt_len);
                else
                    seg_warp = zeros(n_muscles, tgt_len);
                end

                % Place into trial buffer
                if p == 1
                    col_s = 1;
                else
                    col_s = col_e + 1;
                end
                col_e = col_s + tgt_len - 1;
                emg_trial(:, col_s:col_e) = seg_warp;

                % Advance LFP end index (keep tracking in LFP space)
                end_idx_lfp = start_idx_lfp + phase_len_lfp_nom - 1;

            end  % phase loop

            % Append to pool
            pool_emg(:, :, end+1) = emg_trial;

        end  % trial loop

    end  % block loop

    if isempty(pool_emg) || size(pool_emg, 3) == 0
        warning('Subject %s: empty trial pool — skipping.', subj);
        continue;
    end

    % Subject median across trials (mirrors TF v2: median not mean)
    for m = 1:n_muscles
        subj_median_emg{m, s} = median(squeeze(pool_emg(m, :, :)), 2)';  % [1 × tot_length]
    end
    fprintf('  Pooled %d trials.\n', size(pool_emg, 3));

end  % subject loop

%% ===== GRAND AVERAGE =====
ga_emg  = zeros(n_muscles, tot_length);
sem_emg = zeros(n_muscles, tot_length);

for m = 1:n_muscles
    valid = ~cellfun(@isempty, subj_median_emg(m,:));
    if any(valid)
        mat = cell2mat(subj_median_emg(m, valid)');   % [n_valid × tot_length]
        ga_emg(m,:)  = mean(mat, 1);
        sem_emg(m,:) = std(mat, 0, 1) / sqrt(sum(valid));
    end
end

%% ===== FIGURE 1: Grand Average =====
shade = @(ax,t,mn,se,col) fill(ax, [t, fliplr(t)], [mn+se, fliplr(mn-se)], ...
    col, 'EdgeColor','none','FaceAlpha',COL_SHADE,'HandleVisibility','off');

figure('Name','Grand Average EMG Activation (warped)','Position',[50 50 1400 500],'Visible','off');
hold on; grid on;

for m = 1:n_muscles
    shade(gca, time_axis_pct, ga_emg(m,:), sem_emg(m,:), MUSCLE_COLS(m,:));
    plot(time_axis_pct, ga_emg(m,:), 'Color', MUSCLE_COLS(m,:), ...
        'LineWidth', 2.5, 'DisplayName', MUSCLE_LABELS{m});
end

% Phase boundary lines
xline(0,           '--w', 'Mov onset (A)',   'LineWidth', 2, 'LabelVerticalAlignment','bottom');
xline(t_grasp_pct, '--w', 'Grasp start (C)', 'LineWidth', 2, 'LabelVerticalAlignment','bottom');
xline(t_pull_pct,  '--w', 'Pull onset (D)',  'LineWidth', 2, 'LabelVerticalAlignment','bottom');
xline(100,         '--w', 'Mov offset (F)',  'LineWidth', 2, 'LabelVerticalAlignment','bottom');

% Phase background shading
y_lim = [0 max(ga_emg(:) + sem_emg(:)) * 1.15 + 0.01];
for p = 1:n_phases
    if p == 1
        x0 = time_axis_pct(1);
    else
        x0 = time_axis_pct(sum(global_med_len(1:p-1)) + 1);
    end
    x1 = time_axis_pct(sum(global_med_len(1:p)));
    patch([x0 x1 x1 x0], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], PHASE_COLS{p}, ...
        'FaceAlpha', 0.07, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xlabel('Normalised trial time (% of movement)', 'FontSize', 12);
ylabel('Norm. EMG envelope (re 95th pct)', 'FontSize', 12);
ylim(y_lim);
xlim([time_axis_pct(1), time_axis_pct(end)]);
legend('Location','best','FontSize',11);
set(gca,'FontSize',11);
title(sprintf('Grand Average EMG Activation — Warped  (n=%d)', n_subj), 'FontSize', 13);

set(gcf,'Units','pixels','Position',[50 50 1400 500]);
print(fullfile(RESULTS_DIR,'EMG_grandAvg_warped'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'EMG_grandAvg_warped.fig'));

%% ===== FIGURE 2: Per-Subject (one subplot per subject per muscle) =====
figure('Name','Per-Subject EMG Activation (warped)','Position',[50 50 1800 250*n_muscles],'Visible','off');

for m = 1:n_muscles
    for s = 1:n_subj
        ax = subplot(n_muscles, n_subj, (m-1)*n_subj + s);
        hold(ax,'on'); grid(ax,'on');

        if isempty(subj_median_emg{m,s})
            title(ax, [SUBJECTS{s} ' — N/A'], 'FontSize', 8); continue;
        end

        % Phase shading
        ax_ylim = [0 max(subj_median_emg{m,s}) * 1.3 + 0.01];
        for p = 1:n_phases
            if p == 1, x0 = time_axis_pct(1);
            else,       x0 = time_axis_pct(sum(global_med_len(1:p-1))+1); end
            x1 = time_axis_pct(sum(global_med_len(1:p)));
            patch(ax, [x0 x1 x1 x0], [ax_ylim(1) ax_ylim(1) ax_ylim(2) ax_ylim(2)], ...
                PHASE_COLS{p}, 'FaceAlpha',0.1,'EdgeColor','none');
        end

        plot(ax, time_axis_pct, subj_median_emg{m,s}, 'Color', MUSCLE_COLS(m,:), 'LineWidth', 1.5);
        xline(ax, 0,           '--k', 'LineWidth', 1);
        xline(ax, t_grasp_pct, '--k', 'LineWidth', 1);
        xline(ax, t_pull_pct,  '--k', 'LineWidth', 1);
        xline(ax, 100,         '--k', 'LineWidth', 1);

        xlim(ax, [time_axis_pct(1), time_axis_pct(end)]);
        ylim(ax, ax_ylim);
        if m == 1, title(ax, SUBJECTS{s}, 'FontSize', 9); end
        if s == 1, ylabel(ax, MUSCLE_LABELS{m}, 'FontSize', 9); end
        if m == n_muscles, xlabel(ax, '% movement', 'FontSize', 8); end
        set(ax, 'FontSize', 8);
    end
end

sgtitle('Per-Subject Warped EMG Activation', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1800 250*n_muscles]);
print(fullfile(RESULTS_DIR,'EMG_perSubject_warped'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'EMG_perSubject_warped.fig'));
close(gcf);

%% ===== SAVE =====
EMG_WARP.subjects       = SUBJECTS;
EMG_WARP.muscles        = MUSCLE_LABELS;
EMG_WARP.time_axis_pct  = time_axis_pct;
EMG_WARP.global_med_len = global_med_len;
EMG_WARP.fs_lfp         = fs_lfp;
EMG_WARP.ga_emg         = ga_emg;    % [n_muscles × tot_length]
EMG_WARP.sem_emg        = sem_emg;
EMG_WARP.per_subj       = subj_median_emg;   % {n_muscles × n_subj}
EMG_WARP.t_phase_pct    = [0, t_grasp_pct, t_pull_pct, 100];
EMG_WARP.params.bp_band    = [BP_LOW, BP_HIGH];
EMG_WARP.params.env_cutoff = ENV_CUTOFF;
EMG_WARP.params.norm_pctile = NORM_PCTILE;
EMG_WARP.params.rest_offset_s = REST_OFFSET_S;

save(fullfile(RESULTS_DIR,'EMG_activation_warped.mat'), 'EMG_WARP', '-v7.3');
fprintf('\nSaved EMG warped activation to %s\n', RESULTS_DIR);
