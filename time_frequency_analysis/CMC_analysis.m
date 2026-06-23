%% CMC_analysis.m
% =========================================================================
% Corticomuscular Coherence (CMC) between STN LFP and contralateral EMG.
%
% Method: Multi-trial coherence (Halliday et al. 1995 / Rosenberg et al. 1989)
%   For each subject and epoch type:
%     1. Load LFP trials (from Epochs_allSubjects.mat).
%     2. Load continuous EMG, apply BP 20–180 Hz, epoch to same windows.
%     3. For each trial, Hann-taper and FFT both LFP and EMG signals.
%     4. Accumulate cross-spectrum S_xy and auto-spectra S_xx, S_yy.
%     5. Coherence = |sum(S_xy)|^2 / (sum(S_xx) * sum(S_yy))
%        This is unbiased: averaging complex spectra preserves phase
%        cancellation for non-coherent activity.
%     6. 95% confidence bound: 1 - alpha^(1/(N-1)), N = number of trials.
%
% Epoch time windows analysed:
%   Movement epoch:  [0, +1.0] s from onset epoch  (movement period only)
%   Pull epoch:      [0, +0.875] s from pull epoch  (pull period only)
%   Baseline:        [-0.5, 0] s from onset epoch   (pre-movement rest)
%
% Frequency resolution determined by window length:
%   Movement: 1.0 s  → Δf = 1.0 Hz
%   Pull:     0.875 s → Δf ≈ 1.14 Hz
%   Baseline: 0.5 s  → Δf = 2.0 Hz
%
% Outputs (saved to RESULTS_DIR/CMC/):
%   CMC_grandAvg.png / .fig          — coherence spectra per muscle
%   CMC_perSubject_MovOnset.png / .fig
%   CMC_perSubject_Pull.png / .fig
%   CMC_byBand.png / .fig            — bar plot per frequency band
%   CMC_allSubjects.mat
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
LFP_MAT     = 'H:\Parkinson_ReachGrasp\Reprocessing\RESULTS_final\Epochs\Epochs_allSubjects.mat';
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
INPUT_SUB   = fullfile('01_Extracted', 'EMG_KIN');
LFP_SUB     = fullfile('Preprocessed', 'LFP');   % event source (*_wEv.set, at 400 Hz)
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'CMC');

% Epoch windows to analyse (seconds, relative to event)
WIN_MOVE     = [0,    1.000];   % movement phase in onset epoch
WIN_PULL     = [0,    0.875];   % pull phase in pull epoch
WIN_BASELINE = [-0.5, 0.000];  % rest baseline (from onset epoch)

% Frequency axis
FREQ_MIN   = 1;     % Hz — lowest frequency to display
FREQ_MAX   = 100;   % Hz — highest frequency to display

% Contralateral EMG channel selection
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};
TARGET_MUSCLES    = {'iod', 'triceps', 'deltoid'};
MUSCLE_LABELS     = {'IOD', 'Triceps', 'Deltoid'};
n_muscles         = numel(TARGET_MUSCLES);

% EMG preprocessing (same as BatchEMGProcessing_v2)
BP_LOW  = 20;  BP_HIGH = 180;  BP_ORDER = 4;
NOTCH_FREQS = [50, 100];  NOTCH_BW = 2;  NOTCH_ORDER = 2;

% Significance threshold
ALPHA = 0.05;

% Frequency bands for summary bar plot
BANDS = struct( ...
    'name', {'δ (1–4)','θ (4–8)','α (8–13)','β (13–30)','lγ (30–60)'}, ...
    'lim',  {[1 4],    [4 8],    [8 13],    [13 30],    [30 60]});
n_bands = numel(BANDS);

% Display colours per muscle
MUSCLE_COLS = [0.15 0.65 0.25;
               0.80 0.30 0.10;
               0.20 0.40 0.80];

%% ===== BIOSIG STUB FIX =====
if ~isempty(which('eeglab'))
    biosig_stubs = fullfile(fileparts(which('eeglab')), ...
        'plugins', 'Biosig3.8.4', 'biosig', 'maybe-missing');
    if exist(biosig_stubs, 'dir'), rmpath(biosig_stubs); end
end

%% ===== LOAD LFP EPOCHS =====
fprintf('Loading LFP epochs...\n');
load(LFP_MAT, 'EPOCHS');
SUBJECTS = EPOCHS.params.subjects;
n_subj   = numel(SUBJECTS);

eeglab nogui;
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

%% ===== PRE-ALLOCATE =====
% Will hold per-subject coherence [n_freq × n_muscles] for each window
% Stored as structs after first subject sets frequency vector.
coh_move_subj = cell(n_subj, n_muscles);   % coherence during movement
coh_pull_subj = cell(n_subj, n_muscles);
coh_base_subj = cell(n_subj, n_muscles);   % coherence during baseline
ci95_subj     = zeros(n_subj, 3);          % 95% CI bounds: [move, pull, base]

fr_move = [];  fr_pull = [];  fr_base = [];

%% ===== SUBJECT LOOP =====
for s = 1:n_subj
    subj = SUBJECTS{s};
    if ~isfield(EPOCHS, subj)
        warning('Subject %s not in EPOCHS — skipping.', subj);
        continue;
    end

    emg_path  = fullfile(BASE_PATH, subj, INPUT_SUB);
    set_files = dir(fullfile(emg_path, '*.set'));
    if isempty(set_files)
        warning('No EMG .set files for %s — skipping.', subj);
        continue;
    end

    fprintf('\n===== Subject: %s =====\n', subj);
    use_dominant = ismember(subj, DOMINANT_SUBJECTS);

    % LFP *_wEv.set files — used only for their event structure
    lfp_path  = fullfile(BASE_PATH, subj, LFP_SUB);
    lfp_files = dir(fullfile(lfp_path, '*_wEv.set'));
    if isempty(lfp_files)
        warning('No LFP *_wEv.set for %s — cannot source events, skipping.', subj);
        continue;
    end

    % Retrieve LFP epoch data
    t_onset = EPOCHS.(subj).onset.times;   % [1 × 600]
    t_pull  = EPOCHS.(subj).pull.times;    % [1 × 470]
    fs      = EPOCHS.(subj).fs;
    lfp_on  = EPOCHS.(subj).onset.data;    % [600 × n_trials_on]
    lfp_pu  = EPOCHS.(subj).pull.data;     % [470 × n_trials_pu]
    n_t_on  = EPOCHS.(subj).onset.n;
    n_t_pu  = EPOCHS.(subj).pull.n;

    % ---------------------------------------------------------------
    % Epoch the BP-filtered EMG to same windows, across all blocks
    % ---------------------------------------------------------------
    emg_on_ep  = cell(n_muscles, 1);   % each: [n_samp_on × n_trials]
    emg_pu_ep  = cell(n_muscles, 1);

    for f = 1:numel(set_files)
        if f > numel(lfp_files)
            warning('  Block %d: no matching LFP event file — skipping.', f);
            continue;
        end

        EEG    = pop_loadset('filename', set_files(f).name, 'filepath', emg_path);
        fs_emg = EEG.srate;

        % Resample EMG to LFP rate if needed — coherence requires shared fs.
        % Do this BEFORE injecting events so latencies stay in LFP sample space.
        if fs_emg ~= fs
            fprintf('  Block %d: resampling EMG %d→%d Hz\n', f, fs_emg, fs);
            EEG    = pop_resample(EEG, fs);
            fs_emg = fs;
        end

        % ---- Inject kinematic events from preprocessed LFP ----
        % EMG_KIN .set files carry no A_T*/D_T* markers; copy from LFP.
        % After resampling EMG to LFP fs, latencies are directly compatible.
        LFP_ev   = pop_loadset('filename', lfp_files(f).name, 'filepath', lfp_path);
        kin_mask = ~cellfun(@isempty, regexp({LFP_ev.event.type}, '^[A-F]_T\d+$'));
        kin_evs  = LFP_ev.event(kin_mask);
        if isempty(kin_evs)
            warning('  Block %d: LFP file has no kinematic events — skipping.', f);
            continue;
        end
        EEG.event = [EEG.event, kin_evs];
        EEG        = eeg_checkset(EEG);

        % Channel selection
        ch_labels = lower({EEG.chanlocs.labels});
        if use_dominant
            sel_mask = contains(ch_labels, '_dominant_emg');
        else
            sel_mask = contains(ch_labels, '_emg');
        end
        ch_all = find(sel_mask);
        if isempty(ch_all), continue; end

        muscle_ch = zeros(1, n_muscles);
        for m = 1:n_muscles
            hits = find(contains(ch_labels(ch_all), TARGET_MUSCLES{m}));
            if ~isempty(hits), muscle_ch(m) = ch_all(hits(1)); end
        end
        if ~any(muscle_ch), continue; end

        % Build filters
        [b_bp, a_bp] = butter(BP_ORDER, [BP_LOW, BP_HIGH]/(fs/2), 'bandpass');
        notch_filt = cell(numel(NOTCH_FREQS), 2);
        for ni = 1:numel(NOTCH_FREQS)
            fn = NOTCH_FREQS(ni);
            if fn < fs/2
                bw = NOTCH_BW/(fs/2); fc = fn/(fs/2);
                [notch_filt{ni,1}, notch_filt{ni,2}] = ...
                    butter(NOTCH_ORDER, [fc-bw/2, fc+bw/2], 'stop');
            end
        end

        % Preprocess each muscle → store as row in bp_data [n_muscles × n_samp]
        bp_data = zeros(n_muscles, EEG.pnts);
        for m = 1:n_muscles
            if muscle_ch(m) == 0, continue; end
            sig = double(EEG.data(muscle_ch(m), :));
            sig = sig - mean(sig);
            sig = filtfilt(b_bp, a_bp, sig);
            for ni = 1:numel(NOTCH_FREQS)
                if ~isempty(notch_filt{ni,1})
                    sig = filtfilt(notch_filt{ni,1}, notch_filt{ni,2}, sig);
                end
            end
            bp_data(m, :) = sig;
        end

        % Epoch using EEGLAB — reuse event structure from EEG
        EEG_proc        = EEG;
        EEG_proc.data   = bp_data;
        EEG_proc.nbchan = n_muscles;
        for m = 1:n_muscles
            EEG_proc.chanlocs(m).labels = MUSCLE_LABELS{m};
        end
        EEG_proc.chanlocs = EEG_proc.chanlocs(1:n_muscles);

        ev_types = {EEG.event.type};
        types_A  = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^A_T\d+$'))));
        types_D  = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^D_T\d+$'))));

        if ~isempty(types_A)
            EEG_ep = pop_epoch(EEG_proc, types_A, ...
                [-EPOCHS.params.pre_onset_s, EPOCHS.params.post_onset_s]);
            for m = 1:n_muscles
                if size(EEG_ep.data,3) > 0
                    trials = squeeze(EEG_ep.data(m,:,:));
                    if isvector(trials), trials = trials(:); end
                    emg_on_ep{m} = [emg_on_ep{m}, trials];
                end
            end
        end

        if ~isempty(types_D)
            EEG_ep = pop_epoch(EEG_proc, types_D, ...
                [-EPOCHS.params.pre_pull_s, EPOCHS.params.post_pull_s]);
            for m = 1:n_muscles
                if size(EEG_ep.data,3) > 0
                    trials = squeeze(EEG_ep.data(m,:,:));
                    if isvector(trials), trials = trials(:); end
                    emg_pu_ep{m} = [emg_pu_ep{m}, trials];
                end
            end
        end

    end  % block loop

    % Check trial counts
    n_emg_on = size(emg_on_ep{1}, 2);
    n_emg_pu = size(emg_pu_ep{1}, 2);
    if n_emg_on ~= n_t_on
        warning('  Onset trial mismatch: LFP=%d, EMG=%d — using min.', n_t_on, n_emg_on);
    end
    if n_emg_pu ~= n_t_pu
        warning('  Pull trial mismatch: LFP=%d, EMG=%d — using min.', n_t_pu, n_emg_pu);
    end
    n_use_on = min(n_t_on, n_emg_on);
    n_use_pu = min(n_t_pu, n_emg_pu);

    % ---------------------------------------------------------------
    % Extract time windows
    % ---------------------------------------------------------------
    mask_move = t_onset >= WIN_MOVE(1)     & t_onset < WIN_MOVE(2);
    mask_base = t_onset >= WIN_BASELINE(1) & t_onset < WIN_BASELINE(2);
    mask_pull = t_pull  >= WIN_PULL(1)     & t_pull  < WIN_PULL(2);

    n_samp_move = sum(mask_move);
    n_samp_pull = sum(mask_pull);
    n_samp_base = sum(mask_base);

    % ---------------------------------------------------------------
    % Compute multi-trial coherence for each muscle and window
    % ---------------------------------------------------------------
    for m = 1:n_muscles
        emg_on = emg_on_ep{m};   % [n_samp_on × n_trials]
        emg_pu = emg_pu_ep{m};

        % --- Movement window ---
        if n_use_on > 1 && ~isempty(emg_on)
            [coh, fr_m] = multitrial_coherence( ...
                lfp_on(mask_move, 1:n_use_on), ...
                emg_on(mask_move, 1:n_use_on), fs);
            coh_move_subj{s,m} = coh;
            if isempty(fr_move), fr_move = fr_m; end
        end

        % --- Baseline window ---
        if n_use_on > 1 && ~isempty(emg_on)
            [coh, fr_b] = multitrial_coherence( ...
                lfp_on(mask_base, 1:n_use_on), ...
                emg_on(mask_base, 1:n_use_on), fs);
            coh_base_subj{s,m} = coh;
            if isempty(fr_base), fr_base = fr_b; end
        end

        % --- Pull window ---
        if n_use_pu > 1 && ~isempty(emg_pu)
            [coh, fr_p] = multitrial_coherence( ...
                lfp_pu(mask_pull, 1:n_use_pu), ...
                emg_pu(mask_pull, 1:n_use_pu), fs);
            coh_pull_subj{s,m} = coh;
            if isempty(fr_pull), fr_pull = fr_p; end
        end
    end

    % Store 95% CI bound per subject (depends on trial count)
    if n_use_on > 1
        ci95_subj(s,1) = 1 - ALPHA^(1/(n_use_on - 1));
        ci95_subj(s,3) = 1 - ALPHA^(1/(n_use_on - 1));
    end
    if n_use_pu > 1
        ci95_subj(s,2) = 1 - ALPHA^(1/(n_use_pu - 1));
    end

    fprintf('  Done — move/baseline N=%d, pull N=%d\n', n_use_on, n_use_pu);

end  % subject loop

%% ===== GRAND AVERAGE COHERENCE =====
freq_mask_move = fr_move >= FREQ_MIN & fr_move <= FREQ_MAX;
freq_mask_pull = fr_pull >= FREQ_MIN & fr_pull <= FREQ_MAX;
freq_mask_base = fr_base >= FREQ_MIN & fr_base <= FREQ_MAX;

fr_plot_move = fr_move(freq_mask_move);
fr_plot_pull = fr_pull(freq_mask_pull);
fr_plot_base = fr_base(freq_mask_base);

ga_coh_move = zeros(numel(fr_plot_move), n_muscles);
ga_coh_pull = zeros(numel(fr_plot_pull), n_muscles);
ga_coh_base = zeros(numel(fr_plot_base), n_muscles);
sem_coh_move = zeros(numel(fr_plot_move), n_muscles);
sem_coh_pull = zeros(numel(fr_plot_pull), n_muscles);
sem_coh_base = zeros(numel(fr_plot_base), n_muscles);

for m = 1:n_muscles
    valid_m = ~cellfun(@isempty, coh_move_subj(:,m));
    if any(valid_m)
        mat = cell2mat(cellfun(@(c) c(freq_mask_move)', coh_move_subj(valid_m,m), 'UniformOutput',false));
        ga_coh_move(:,m)  = mean(mat, 1)';
        sem_coh_move(:,m) = std(mat, 0, 1)' / sqrt(sum(valid_m));
    end
    valid_p = ~cellfun(@isempty, coh_pull_subj(:,m));
    if any(valid_p)
        mat = cell2mat(cellfun(@(c) c(freq_mask_pull)', coh_pull_subj(valid_p,m), 'UniformOutput',false));
        ga_coh_pull(:,m)  = mean(mat, 1)';
        sem_coh_pull(:,m) = std(mat, 0, 1)' / sqrt(sum(valid_p));
    end
    valid_b = ~cellfun(@isempty, coh_base_subj(:,m));
    if any(valid_b)
        mat = cell2mat(cellfun(@(c) c(freq_mask_base)', coh_base_subj(valid_b,m), 'UniformOutput',false));
        ga_coh_base(:,m)  = mean(mat, 1)';
        sem_coh_base(:,m) = std(mat, 0, 1)' / sqrt(sum(valid_b));
    end
end

% Grand average 95% CI (median across subjects)
ci_move = median(ci95_subj(ci95_subj(:,1)>0, 1));
ci_pull = median(ci95_subj(ci95_subj(:,2)>0, 2));
ci_base = median(ci95_subj(ci95_subj(:,3)>0, 3));

%% ===== FIGURE 1: Grand Average CMC — Movement vs. Baseline =====
figure('Name','Grand Average CMC — Movement vs Baseline', ...
    'Position',[50 50 1500 500*n_muscles],'Visible','off');

n_rows = n_muscles;
shade = @(ax,f,mn,se,col) fill(ax, [f(:)', fliplr(f(:)')], ...
    [mn(:)'+se(:)', fliplr(mn(:)'-se(:)')], ...
    col, 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');

for m = 1:n_muscles
    ax = subplot(n_muscles, 2, (m-1)*2 + 1);
    hold on; grid on;
    shade(ax, fr_plot_base, ga_coh_base(:,m)', sem_coh_base(:,m)', [0.5 0.5 0.5]);
    shade(ax, fr_plot_move, ga_coh_move(:,m)', sem_coh_move(:,m)', MUSCLE_COLS(m,:));
    plot(ax, fr_plot_base, ga_coh_base(:,m), 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.5, 'DisplayName', 'Baseline');
    plot(ax, fr_plot_move, ga_coh_move(:,m), 'Color', MUSCLE_COLS(m,:), ...
        'LineWidth', 2,   'DisplayName', 'Movement');
    yline(ax, ci_move, '--k', '95% CI', 'LineWidth', 1, 'FontSize', 8);
    draw_band_lines(ax, BANDS, FREQ_MIN, FREQ_MAX);
    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Coherence','FontSize',10);
    title(sprintf('%s — Mov vs Baseline', MUSCLE_LABELS{m}),'FontSize',11);
    legend('Location','northeast','FontSize',8);
    set(ax,'FontSize',9); xlim([FREQ_MIN FREQ_MAX]); ylim([0 1]);
    patch(ax, [13 30 30 13],[0 0 1 1],[0.9 0.9 0.2],'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');

    ax = subplot(n_muscles, 2, (m-1)*2 + 2);
    hold on; grid on;
    shade(ax, fr_plot_base, ga_coh_base(:,m)', sem_coh_base(:,m)', [0.5 0.5 0.5]);
    shade(ax, fr_plot_pull, ga_coh_pull(:,m)', sem_coh_pull(:,m)', MUSCLE_COLS(m,:));
    plot(ax, fr_plot_base, ga_coh_base(:,m), 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.5, 'DisplayName', 'Baseline');
    plot(ax, fr_plot_pull, ga_coh_pull(:,m), 'Color', MUSCLE_COLS(m,:), ...
        'LineWidth', 2,   'DisplayName', 'Pull');
    yline(ax, ci_pull, '--k', '95% CI', 'LineWidth', 1, 'FontSize', 8);
    draw_band_lines(ax, BANDS, FREQ_MIN, FREQ_MAX);
    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Coherence','FontSize',10);
    title(sprintf('%s — Pull vs Baseline', MUSCLE_LABELS{m}),'FontSize',11);
    legend('Location','northeast','FontSize',8);
    set(ax,'FontSize',9); xlim([FREQ_MIN FREQ_MAX]); ylim([0 1]);
    patch(ax, [13 30 30 13],[0 0 1 1],[0.9 0.9 0.2],'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
end

sgtitle(sprintf('STN LFP–EMG Coherence — Grand Average  (n=%d)', n_subj), 'FontSize', 14);
set(gcf,'Units','pixels','Position',[50 50 1500 500*n_muscles]);
print(fullfile(RESULTS_DIR,'CMC_grandAvg'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'CMC_grandAvg.fig'));

%% ===== FIGURE 2: Band-averaged CMC bar plot =====
% For each muscle × band: mean CMC during movement vs baseline
figure('Name','CMC by Frequency Band','Position',[50 50 1400 400],'Visible','off');

for m = 1:n_muscles
    ax = subplot(1, n_muscles, m);
    hold on; grid on;

    band_move = zeros(1, n_bands);
    band_base = zeros(1, n_bands);
    for b = 1:n_bands
        fmask_m = fr_plot_move >= BANDS(b).lim(1) & fr_plot_move <= BANDS(b).lim(2);
        fmask_b = fr_plot_base >= BANDS(b).lim(1) & fr_plot_base <= BANDS(b).lim(2);
        if any(fmask_m), band_move(b) = mean(ga_coh_move(fmask_m, m)); end
        if any(fmask_b), band_base(b) = mean(ga_coh_base(fmask_b, m)); end
    end

    x = 1:n_bands;
    bar(ax, x - 0.2, band_base, 0.35, 'FaceColor', [0.6 0.6 0.6], 'DisplayName', 'Baseline');
    bar(ax, x + 0.2, band_move, 0.35, 'FaceColor', MUSCLE_COLS(m,:),    'DisplayName', 'Movement');

    set(ax, 'XTick', x, 'XTickLabel', {BANDS.name}, 'FontSize', 9);
    xtickangle(30);
    ylabel('Mean Coherence','FontSize',10);
    title(MUSCLE_LABELS{m},'FontSize',11);
    if m == 1, legend('Location','northwest','FontSize',8); end
    ylim([0 max([band_move, band_base]) * 1.3 + 0.01]);
end

sgtitle('STN LFP–EMG Coherence by Frequency Band (grand avg ± SEM)', 'FontSize', 13);
set(gcf,'Units','pixels','Position',[50 50 1400 400]);
print(fullfile(RESULTS_DIR,'CMC_byBand'), '-dpng', '-r300');
savefig(fullfile(RESULTS_DIR,'CMC_byBand.fig'));

%% ===== FIGURES 3–4: Per-Subject CMC =====
for cond = 1:2
    if cond == 1
        coh_cell = coh_move_subj;
        fr_plot  = fr_plot_move;
        ci_line  = ci_move;
        cond_lbl = 'Movement';
        fig_name = 'CMC_perSubject_MovOnset';
    else
        coh_cell = coh_pull_subj;
        fr_plot  = fr_plot_pull;
        ci_line  = ci_pull;
        cond_lbl = 'Pull';
        fig_name = 'CMC_perSubject_Pull';
    end

    figure('Name', sprintf('Per-Subject CMC — %s', cond_lbl), ...
        'Position',[50 50 1600 300*n_muscles],'Visible','off');

    for m = 1:n_muscles
        for s = 1:n_subj
            ax = subplot(n_muscles, n_subj, (m-1)*n_subj + s);
            hold on; grid on;
            subj = SUBJECTS{s};

            coh_s = coh_cell{s,m};
            if isempty(coh_s)
                title([subj ' — MISSING'],'FontSize',7); continue;
            end

            plot(ax, fr_plot, coh_s(freq_mask_move | freq_mask_pull | freq_mask_base), ...
                'Color', MUSCLE_COLS(m,:), 'LineWidth', 1);

            % use correct mask per condition
            if cond == 1
                plot(ax, fr_plot, coh_s(freq_mask_move), 'Color', MUSCLE_COLS(m,:), 'LineWidth', 1);
            else
                plot(ax, fr_plot, coh_s(freq_mask_pull), 'Color', MUSCLE_COLS(m,:), 'LineWidth', 1);
            end

            yline(ax, ci_line, '--k', 'LineWidth', 0.8);
            set(ax,'FontSize',7); xlim([FREQ_MIN FREQ_MAX]); ylim([0 0.8]);
            if m == 1, title(subj,'FontSize',8); end
            if s == 1, ylabel(MUSCLE_LABELS{m},'FontSize',8); end
            if m == n_muscles, xlabel('Hz','FontSize',7); end
        end
    end

    sgtitle(sprintf('Per-Subject STN LFP–EMG Coherence — %s', cond_lbl), 'FontSize', 13);
    set(gcf,'Units','pixels','Position',[50 50 1600 300*n_muscles]);
    print(fullfile(RESULTS_DIR, fig_name), '-dpng', '-r300');
    savefig(fullfile(RESULTS_DIR, [fig_name '.fig']));
    close(gcf);
end

%% ===== SAVE =====
CMC.subjects         = SUBJECTS;
CMC.muscles          = MUSCLE_LABELS;
CMC.fr_move          = fr_move;
CMC.fr_pull          = fr_pull;
CMC.fr_base          = fr_base;
CMC.ga_coh_move      = ga_coh_move;    % [n_freq × n_muscles]
CMC.ga_coh_pull      = ga_coh_pull;
CMC.ga_coh_base      = ga_coh_base;
CMC.sem_coh_move     = sem_coh_move;
CMC.sem_coh_pull     = sem_coh_pull;
CMC.sem_coh_base     = sem_coh_base;
CMC.per_subj_move    = coh_move_subj;  % {n_subj × n_muscles} cell
CMC.per_subj_pull    = coh_pull_subj;
CMC.per_subj_base    = coh_base_subj;
CMC.ci95_subj        = ci95_subj;
CMC.ci_move          = ci_move;
CMC.ci_pull          = ci_pull;
CMC.ci_base          = ci_base;
CMC.params.win_move     = WIN_MOVE;
CMC.params.win_pull     = WIN_PULL;
CMC.params.win_baseline = WIN_BASELINE;
CMC.params.bp_band      = [BP_LOW, BP_HIGH];
CMC.params.alpha        = ALPHA;

save(fullfile(RESULTS_DIR, 'CMC_allSubjects.mat'), 'CMC', '-v7.3');
fprintf('\nSaved CMC results to %s\n', RESULTS_DIR);

%% ===== LOCAL FUNCTIONS =====

function [coh, fr] = multitrial_coherence(lfp_mat, emg_mat, fs)
% Halliday 1995 multi-trial coherence.
% lfp_mat, emg_mat: [n_samp × n_trials]
% coh: [n_freq × 1] magnitude-squared coherence
% fr:  [n_freq × 1] frequency vector (Hz)

[n_samp, n_trials] = size(lfp_mat);
nfft  = n_samp;                  % use full window length → max freq resolution
win   = hann(n_samp)';           % Hann taper to reduce spectral leakage

Sxx = zeros(nfft/2+1, 1);
Syy = zeros(nfft/2+1, 1);
Sxy = zeros(nfft/2+1, 1);       % complex cross-spectrum

for t = 1:n_trials
    x = lfp_mat(:,t)' .* win;
    y = emg_mat(:,t)' .* win;
    X = fft(x, nfft);
    Y = fft(y, nfft);
    X = X(1:nfft/2+1);
    Y = Y(1:nfft/2+1);
    Sxx = Sxx + real(X .* conj(X));
    Syy = Syy + real(Y .* conj(Y));
    Sxy = Sxy + X .* conj(Y);       % accumulate complex cross-spectrum
end

% Halliday coherence: use averaged complex cross-spectrum
coh = abs(Sxy).^2 ./ (Sxx .* Syy);
coh = coh(:);

fr  = (0 : nfft/2)' * (fs / nfft);
end


function draw_band_lines(ax, bands, fmin, fmax)
% Draw thin vertical dotted lines at canonical band boundaries.
cols = {[0.5 0.3 0.1],[0.5 0.0 0.7],[0.0 0.5 0.0],[0.8 0.0 0.0],[0.6 0.0 0.0]};
limits = unique([bands.lim]);
for k = 1:numel(limits)
    if limits(k) > fmin && limits(k) < fmax
        xline(ax, limits(k), ':', 'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 0.8);
    end
end
end
