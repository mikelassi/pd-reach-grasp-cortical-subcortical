%% CMC_EEG_analysis.m
% =========================================================================
% Corticomuscular Coherence (CMC) between scalp EEG and contralateral EMG.
%
% Method: Multi-trial coherence (Halliday et al. 1995)
%   Same estimator as CMC_analysis.m (STN LFP version) — see that file
%   for full methodological notes. Here applied across all EEG channels,
%   giving a spatial (topoplot) map of scalp CMC per frequency band.
%
% Data sources:
%   EEG  — Preprocessed/EEG/*_manual.set
%          Manually cleaned, kinematic events (A_T1 … F_Tn) embedded.
%          Channel locations assigned from standard_1005.elc via DIPFIT.
%   EMG  — 01_Extracted/EMG_KIN/*.set
%          Contralateral muscles: IOD, Triceps, Deltoid.
%          Band-pass 20–180 Hz applied inline; NOT enveloped.
%
% Note on sampling rates:
%   EEG and EMG may differ (EEG is often at 512/256 Hz, EMG at 2 kHz).
%   The script resamples the higher-rate signal to match the lower rate
%   before computing coherence to keep a common frequency axis.
%
% Epoch windows analysed:
%   Movement: [0, +1.0] s  from movement-onset epoch
%   Pull:     [0, +0.875] s from pull-onset epoch
%   Baseline: [-0.5, 0] s  from movement-onset epoch
%
% Outputs (RESULTS_DIR/CMC_EEG/):
%   CMC_EEG_topo_movement.png / .fig   — topoplots per muscle × band (movement)
%   CMC_EEG_topo_pull.png / .fig       — topoplots per muscle × band (pull)
%   CMC_EEG_topo_difference.png / .fig — topoplots per muscle × band (move−base)
%   CMC_EEG_beta_perSubject.png        — per-subject beta-band topoplots
%   CMC_EEG_spectrum_peakElec.png      — coherence spectrum at peak electrode
%   CMC_EEG_allSubjects.mat
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
SUBJECTS    = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
EEG_SUB     = fullfile('Preprocessed', 'EEG');    % *_manual.set with events
EMG_SUB     = fullfile('01_Extracted', 'EMG_KIN');
LFP_SUB     = fullfile('Preprocessed', 'LFP');   % event source for EMG (*_wEv.set, 400 Hz)
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'CMC_EEG');

% Epoch windows (same as CMC_analysis.m)
WIN_MOVE     = [0,    1.000];
WIN_PULL     = [0,    0.875];
WIN_BASELINE = [-0.5, 0.000];

% EEG epoch windows (from movement onset and pull onset)
PRE_ONSET_S  = 0.5;   POST_ONSET_S = 1.0;
PRE_PULL_S   = 0.3;   POST_PULL_S  = 0.875;

% EMG preprocessing
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};
TARGET_MUSCLES    = {'iod', 'triceps', 'deltoid'};
MUSCLE_LABELS     = {'IOD', 'Triceps', 'Deltoid'};
n_muscles         = numel(TARGET_MUSCLES);
BP_LOW  = 20;  BP_HIGH = 180;  BP_ORDER = 4;
NOTCH_FREQS = [50, 100];  NOTCH_BW = 2;  NOTCH_ORDER = 2;

% Frequency display limits for spectra
FREQ_MIN = 1;
FREQ_MAX = 100;

% Significance threshold
ALPHA = 0.05;

% Frequency bands
BANDS = struct( ...
    'name', {'δ (1–4)','θ (4–8)','α (8–13)','β (13–30)','lγ (30–60)'}, ...
    'lim',  {[1 4],    [4 8],    [8 13],    [13 30],    [30 60]});
n_bands = numel(BANDS);

% wue06 block 3 EEG is corrupted — skip it
SKIP_BLOCKS = containers.Map({'wue06'}, {3});

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

%% ===== FIRST PASS: discover common EEG channels across all subjects =====
fprintf('First pass: discovering EEG channel sets...\n');
ch_sets = cell(1, n_subj);
chanlocs_ref = [];

for s = 1:n_subj
    subj    = SUBJECTS{s};
    eeg_path = fullfile(BASE_PATH, subj, EEG_SUB);
    files    = dir(fullfile(eeg_path, '*_manual.set'));
    if isempty(files), continue; end
    EEG = pop_loadset('filename', files(1).name, 'filepath', eeg_path);
    ch_sets{s} = lower({EEG.chanlocs.labels});
    if isempty(chanlocs_ref)
        chanlocs_ref = EEG.chanlocs;
        fs_eeg_ref   = EEG.srate;
    end
    fprintf('  %s: %d channels, fs=%d Hz\n', subj, EEG.nbchan, EEG.srate);
end

% Intersection of channel labels across all subjects
valid_sets = ch_sets(~cellfun(@isempty, ch_sets));
common_labels = valid_sets{1};
for k = 2:numel(valid_sets)
    common_labels = intersect(common_labels, valid_sets{k});
end
fprintf('\nCommon channels across all subjects: %d\n', numel(common_labels));

% Build common chanlocs from first-subject reference
common_idx_ref = find(ismember(lower({chanlocs_ref.labels}), common_labels));
common_chanlocs = chanlocs_ref(common_idx_ref);
n_ch_common     = numel(common_chanlocs);

%% ===== PRE-ALLOCATE RESULTS =====
% Per-subject coherence: {n_subj × n_muscles}, each [n_ch_common × n_freq]
coh_move_subj = cell(n_subj, n_muscles);
coh_pull_subj = cell(n_subj, n_muscles);
coh_base_subj = cell(n_subj, n_muscles);
ci95_subj     = zeros(n_subj, 3);   % [move, pull, base]

fr_move = [];  fr_pull = [];  fr_base = [];

%% ===== SUBJECT LOOP =====
for s = 1:n_subj
    subj     = SUBJECTS{s};
    eeg_path = fullfile(BASE_PATH, subj, EEG_SUB);
    emg_path = fullfile(BASE_PATH, subj, EMG_SUB);
    eeg_files = dir(fullfile(eeg_path, '*_manual.set'));
    emg_files = dir(fullfile(emg_path, '*.set'));

    if isempty(eeg_files) || isempty(emg_files)
        warning('Missing EEG or EMG files for %s — skipping.', subj);
        continue;
    end

    % Remove corrupted blocks
    if isKey(SKIP_BLOCKS, subj)
        skip_b = SKIP_BLOCKS(subj);
        keep = setdiff(1:numel(eeg_files), skip_b);
        eeg_files = eeg_files(keep);
        emg_files = emg_files(keep);
        fprintf('  %s: skipping block %d (corrupted)\n', subj, skip_b);
    end

    fprintf('\n===== Subject: %s  (%d EEG blocks, %d EMG blocks) =====\n', ...
        subj, numel(eeg_files), numel(emg_files));

    use_dominant = ismember(subj, DOMINANT_SUBJECTS);

    % LFP files as event source — EMG_KIN .set files have no kinematic markers
    lfp_path  = fullfile(BASE_PATH, subj, LFP_SUB);
    lfp_files = dir(fullfile(lfp_path, '*_wEv.set'));
    if isempty(lfp_files)
        warning('No LFP *_wEv.set for %s — cannot source EMG events, skipping.', subj);
        continue;
    end

    % Epoch accumulators: each {m} = [n_ch × n_samp × n_trials]
    eeg_on_ep   = cell(1, 1);        % just one cell, 3D array inside
    eeg_pu_ep   = cell(1, 1);
    emg_on_ep   = cell(n_muscles, 1);
    emg_pu_ep   = cell(n_muscles, 1);
    eeg_on_3d   = [];  eeg_pu_3d   = [];

    fs_eeg = NaN;  fs_emg = NaN;

    for f = 1:min(numel(eeg_files), numel(emg_files))

        % ---- Load EEG ----
        EEG = pop_loadset('filename', eeg_files(f).name, 'filepath', eeg_path);
        if isnan(fs_eeg), fs_eeg = EEG.srate; end

        % Select common channels in this subject's file
        subj_labels = lower({EEG.chanlocs.labels});
        [~, subj_idx] = ismember(common_labels, subj_labels);
        missing = common_labels(subj_idx == 0);
        if ~isempty(missing)
            warning('  Block %d: %d common channels missing in EEG — using available.', ...
                f, numel(missing));
        end
        subj_idx = subj_idx(subj_idx > 0);
        EEG_sel = pop_select(EEG, 'channel', subj_idx);

        % ---- Load and preprocess EMG ----
        EMG = pop_loadset('filename', emg_files(f).name, 'filepath', emg_path);
        if isnan(fs_emg), fs_emg = EMG.srate; end

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
            warning('  Block %d: no target EMG channels for %s — skipping.', f, subj);
            continue;
        end

        % BP filter EMG
        [b_bp, a_bp] = butter(BP_ORDER, [BP_LOW, BP_HIGH]/(fs_emg/2), 'bandpass');
        notch_b = cell(numel(NOTCH_FREQS),2);
        for ni = 1:numel(NOTCH_FREQS)
            fn = NOTCH_FREQS(ni);
            if fn < fs_emg/2
                bw = NOTCH_BW/(fs_emg/2); fc = fn/(fs_emg/2);
                [notch_b{ni,1}, notch_b{ni,2}] = butter(NOTCH_ORDER, [fc-bw/2,fc+bw/2], 'stop');
            end
        end

        bp_emg = zeros(n_muscles, EMG.pnts);
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
            bp_emg(m,:) = sig;
        end

        % ---- Resample if EEG/EMG sampling rates differ ----
        % Downsample both to the lower rate to keep coherence up to min(fs/2)
        fs_coh = min(fs_eeg, fs_emg);
        if fs_eeg ~= fs_coh
            EEG_sel = pop_resample(EEG_sel, fs_coh);
        end
        if fs_emg ~= fs_coh
            bp_emg_ds = resample(bp_emg', round(EMG.pnts * fs_coh / fs_emg), EMG.pnts)';
        else
            bp_emg_ds = bp_emg;
        end

        % Put resampled EMG into an EEG struct for pop_epoch
        EMG_proc          = EMG;
        EMG_proc.data     = bp_emg_ds;
        EMG_proc.nbchan   = n_muscles;
        EMG_proc.pnts     = size(bp_emg_ds, 2);
        EMG_proc.srate    = fs_coh;
        EMG_proc.times    = (0:EMG_proc.pnts-1) / fs_coh * 1000;
        for m = 1:n_muscles
            EMG_proc.chanlocs(m).labels = MUSCLE_LABELS{m};
        end
        EMG_proc.chanlocs = EMG_proc.chanlocs(1:n_muscles);

        % ---- Inject kinematic events from preprocessed LFP into EMG_proc ----
        % EMG_KIN .set files have no A_T*/D_T* markers; copy from LFP.
        % LFP events are at 400 Hz; EMG_proc is at fs_coh — convert.
        if f <= numel(lfp_files)
            LFP_ev   = pop_loadset('filename', lfp_files(f).name, 'filepath', lfp_path);
            kin_mask = ~cellfun(@isempty, regexp({LFP_ev.event.type}, '^[A-F]_T\d+$'));
            kin_evs  = LFP_ev.event(kin_mask);
            fs_lfp   = LFP_ev.srate;
            if ~isempty(kin_evs)
                for k = 1:numel(kin_evs)
                    kin_evs(k).latency = round(kin_evs(k).latency * fs_coh / fs_lfp);
                end
                EMG_proc.event = [EMG_proc.event, kin_evs];
                EMG_proc        = eeg_checkset(EMG_proc);
            end
        end

        % ---- Epoch ----
        % EEG: use events already embedded in the *_manual.set file (added by A03).
        ev_types = {EEG.event.type};
        types_A = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^A_T\d+$'))));
        types_D = unique(ev_types(~cellfun(@isempty, regexp(ev_types, '^D_T\d+$'))));

        if ~isempty(types_A)
            EEG_ep = pop_epoch(EEG_sel, types_A, [-PRE_ONSET_S, POST_ONSET_S]);
            EMG_ep = pop_epoch(EMG_proc, types_A, [-PRE_ONSET_S, POST_ONSET_S]);

            n_use = min(EEG_ep.trials, EMG_ep.trials);
            if n_use > 0
                eeg_block = double(EEG_ep.data(:,:,1:n_use));  % [n_ch × n_samp × n_use]
                if isempty(eeg_on_3d)
                    eeg_on_3d = eeg_block;
                else
                    eeg_on_3d = cat(3, eeg_on_3d, eeg_block);
                end
                for m = 1:n_muscles
                    trials = squeeze(double(EMG_ep.data(m,:,1:n_use)));
                    if isvector(trials), trials = trials(:); end
                    emg_on_ep{m} = [emg_on_ep{m}, trials];
                end
            end
        end

        if ~isempty(types_D)
            EEG_ep = pop_epoch(EEG_sel, types_D, [-PRE_PULL_S, POST_PULL_S]);
            if ~isempty(types_D)
                EMG_ep = pop_epoch(EMG_proc, types_D, [-PRE_PULL_S, POST_PULL_S]);
                n_use = min(EEG_ep.trials, EMG_ep.trials);
                if n_use > 0
                    eeg_block = double(EEG_ep.data(:,:,1:n_use));
                    if isempty(eeg_pu_3d)
                        eeg_pu_3d = eeg_block;
                    else
                        eeg_pu_3d = cat(3, eeg_pu_3d, eeg_block);
                    end
                    for m = 1:n_muscles
                        trials = squeeze(double(EMG_ep.data(m,:,1:n_use)));
                        if isvector(trials), trials = trials(:); end
                        emg_pu_ep{m} = [emg_pu_ep{m}, trials];
                    end
                end
            end
        end

    end  % block loop

    if isempty(eeg_on_3d)
        warning('Subject %s: no valid EEG onset epochs — skipping.', subj);
        continue;
    end

    % ---- Extract analysis windows from epochs ----
    t_onset_eeg = (0 : size(eeg_on_3d,2)-1) / fs_coh - PRE_ONSET_S;
    t_pull_eeg  = (0 : size(eeg_pu_3d,2)-1) / fs_coh - PRE_PULL_S;

    mask_move = t_onset_eeg >= WIN_MOVE(1)     & t_onset_eeg < WIN_MOVE(2);
    mask_base = t_onset_eeg >= WIN_BASELINE(1) & t_onset_eeg < WIN_BASELINE(2);
    mask_pull = t_pull_eeg  >= WIN_PULL(1)     & t_pull_eeg  < WIN_PULL(2);

    n_t_on = size(eeg_on_3d, 3);
    n_t_pu = isempty(eeg_pu_3d) + ~isempty(eeg_pu_3d) * size(eeg_pu_3d, 3);

    fprintf('  Epochs: onset=%d, pull=%d\n', n_t_on, n_t_pu);

    % ---- Multichannel multi-trial coherence per muscle ----
    for m = 1:n_muscles
        emg_on = emg_on_ep{m};   % [n_samp_on × n_trials]
        emg_pu = emg_pu_ep{m};

        n_use_on = min(n_t_on, size(emg_on,2));
        if n_use_on > 1
            % Movement window
            [coh_m, fr_m] = mc_multitrial_coh( ...
                eeg_on_3d(:, mask_move, 1:n_use_on), ...
                emg_on(mask_move, 1:n_use_on), fs_coh);
            coh_move_subj{s,m} = coh_m;   % [n_ch_common × n_freq]
            if isempty(fr_move), fr_move = fr_m; end

            % Baseline window
            [coh_b, fr_b] = mc_multitrial_coh( ...
                eeg_on_3d(:, mask_base, 1:n_use_on), ...
                emg_on(mask_base, 1:n_use_on), fs_coh);
            coh_base_subj{s,m} = coh_b;
            if isempty(fr_base), fr_base = fr_b; end

            ci95_subj(s,1) = 1 - ALPHA^(1/(n_use_on - 1));
            ci95_subj(s,3) = ci95_subj(s,1);
        end

        if ~isempty(eeg_pu_3d) && ~isempty(emg_pu)
            n_use_pu = min(size(eeg_pu_3d,3), size(emg_pu,2));
            if n_use_pu > 1
                [coh_p, fr_p] = mc_multitrial_coh( ...
                    eeg_pu_3d(:, mask_pull, 1:n_use_pu), ...
                    emg_pu(mask_pull, 1:n_use_pu), fs_coh);
                coh_pull_subj{s,m} = coh_p;
                if isempty(fr_pull), fr_pull = fr_p; end
                ci95_subj(s,2) = 1 - ALPHA^(1/(n_use_pu - 1));
            end
        end
    end

    fprintf('  Done.\n');
end  % subject loop

%% ===== GRAND AVERAGE COHERENCE =====
fmask_move = fr_move >= FREQ_MIN & fr_move <= FREQ_MAX;
fmask_pull = fr_pull >= FREQ_MIN & fr_pull <= FREQ_MAX;
fmask_base = fr_base >= FREQ_MIN & fr_base <= FREQ_MAX;

fr_plot_m = fr_move(fmask_move);
fr_plot_p = fr_pull(fmask_pull);
fr_plot_b = fr_base(fmask_base);

n_fm = sum(fmask_move);
n_fp = sum(fmask_pull);

ga_move = zeros(n_ch_common, n_fm, n_muscles);
ga_pull = zeros(n_ch_common, n_fp, n_muscles);
ga_base = zeros(n_ch_common, sum(fmask_base), n_muscles);
n_valid_move = zeros(n_muscles, 1);
n_valid_pull = zeros(n_muscles, 1);

for m = 1:n_muscles
    valid_m = find(~cellfun(@isempty, coh_move_subj(:,m)));
    if ~isempty(valid_m)
        stack = cat(3, coh_move_subj{valid_m, m});   % [n_ch × n_freq × n_valid]
        ga_move(:,:,m)   = mean(stack(:, fmask_move, :), 3);
        n_valid_move(m) = numel(valid_m);
    end
    valid_p = find(~cellfun(@isempty, coh_pull_subj(:,m)));
    if ~isempty(valid_p)
        stack = cat(3, coh_pull_subj{valid_p, m});
        ga_pull(:,:,m)   = mean(stack(:, fmask_pull, :), 3);
        n_valid_pull(m) = numel(valid_p);
    end
    valid_b = find(~cellfun(@isempty, coh_base_subj(:,m)));
    if ~isempty(valid_b)
        stack = cat(3, coh_base_subj{valid_b, m});
        ga_base(:,:,m) = mean(stack(:, fmask_base, :), 3);
    end
end

ci_move = median(ci95_subj(ci95_subj(:,1)>0, 1));
ci_pull = median(ci95_subj(ci95_subj(:,2)>0, 2));

% Band-averaged coherence: [n_ch × n_bands × n_muscles] for move/pull/base
band_coh_move = band_average(ga_move, fr_plot_m, BANDS);
band_coh_pull = band_average(ga_pull, fr_plot_p, BANDS);
band_coh_base = band_average(ga_base, fr_plot_b, BANDS);
band_coh_diff = band_coh_move - band_coh_base;   % movement vs rest

%% ===== FIGURE 1: Grand Average Topoplots — Movement =====
fprintf('\nPlotting topoplots...\n');
make_topo_figure(band_coh_move, common_chanlocs, BANDS, MUSCLE_LABELS, ...
    'Movement', 'hot', [0 max(band_coh_move(:))], ...
    fullfile(RESULTS_DIR, 'CMC_EEG_topo_movement'), n_valid_move);

%% ===== FIGURE 2: Grand Average Topoplots — Movement minus Baseline =====
lim_diff = max(abs(band_coh_diff(:)));
make_topo_figure(band_coh_diff, common_chanlocs, BANDS, MUSCLE_LABELS, ...
    'Movement − Baseline', 'rdbu', [-lim_diff lim_diff], ...
    fullfile(RESULTS_DIR, 'CMC_EEG_topo_difference'), n_valid_move);

%% ===== FIGURE 3: Grand Average Topoplots — Pull =====
make_topo_figure(band_coh_pull, common_chanlocs, BANDS, MUSCLE_LABELS, ...
    'Pull', 'hot', [0 max(band_coh_pull(:))], ...
    fullfile(RESULTS_DIR, 'CMC_EEG_topo_pull'), n_valid_pull);

%% ===== FIGURE 4: Per-Subject Beta-Band Topoplots =====
% Find beta band index
beta_band_idx = find(strcmp({BANDS.name}, 'β (13–30)'));
if isempty(beta_band_idx), beta_band_idx = 4; end

figure('Name','Per-Subject Beta CMC','Position',[50 50 1600 300*n_muscles],'Visible','off');
cmax_beta = 0;
for s = 1:n_subj
    for m = 1:n_muscles
        c = coh_move_subj{s,m};
        if ~isempty(c)
            band_vals = mean(c(:, fr_move >= BANDS(beta_band_idx).lim(1) & ...
                               fr_move <= BANDS(beta_band_idx).lim(2)), 2);
            cmax_beta = max(cmax_beta, max(band_vals));
        end
    end
end

for m = 1:n_muscles
    for s = 1:n_subj
        ax = subplot(n_muscles, n_subj, (m-1)*n_subj + s);
        c = coh_move_subj{s,m};
        if isempty(c)
            axis(ax,'off'); text(0.5,0.5,'N/A','HorizontalAlignment','center','Parent',ax);
            continue;
        end
        fmask_beta = fr_move >= BANDS(beta_band_idx).lim(1) & ...
                     fr_move <= BANDS(beta_band_idx).lim(2);
        beta_vals = mean(c(:, fmask_beta), 2);   % [n_ch × 1]
        topoplot(beta_vals, common_chanlocs, ...
            'maplimits', [0 cmax_beta + eps], ...
            'electrodes', 'off', 'colormap', hot(256), 'style', 'map');
        if m == 1, title(SUBJECTS{s}, 'FontSize', 7); end
        if s == 1, ylabel(ax, MUSCLE_LABELS{m}, 'FontSize', 9); end
    end
end

sgtitle(sprintf('Per-Subject Beta CMC (13–30 Hz) — Movement  (EEG vs EMG)'), 'FontSize', 12);
set(gcf,'Units','pixels','Position',[50 50 1600 300*n_muscles]);
colormap(hot(256));
print(fullfile(RESULTS_DIR,'CMC_EEG_beta_perSubject'), '-dpng', '-r200');
savefig(fullfile(RESULTS_DIR,'CMC_EEG_beta_perSubject.fig'));
close(gcf);

%% ===== FIGURE 5: Coherence Spectrum at Peak Electrode =====
figure('Name','Peak-Electrode CMC Spectrum','Position',[50 50 1400 420],'Visible','off');
MUSCLE_COLS = [0.15 0.65 0.25; 0.80 0.30 0.10; 0.20 0.40 0.80];

for m = 1:n_muscles
    ax = subplot(1, n_muscles, m);
    hold on; grid on;

    % Find electrode with highest beta CMC during movement
    fmask_beta = fr_plot_m >= BANDS(beta_band_idx).lim(1) & ...
                 fr_plot_m <= BANDS(beta_band_idx).lim(2);
    beta_avg_per_ch = mean(ga_move(:, fmask_beta, m), 2);
    [~, peak_ch_idx] = max(beta_avg_per_ch);
    peak_label = common_chanlocs(peak_ch_idx).labels;

    plot(ax, fr_plot_b, squeeze(ga_base(peak_ch_idx, :, m)), ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'Baseline');
    plot(ax, fr_plot_m, squeeze(ga_move(peak_ch_idx, :, m)), ...
        'Color', MUSCLE_COLS(m,:), 'LineWidth', 2, 'DisplayName', 'Movement');

    yline(ci_move, '--k', '95% CI', 'LineWidth', 1, 'FontSize', 8);
    patch(ax, [13 30 30 13],[0 0 1 1],[0.9 0.9 0.2], ...
        'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');

    xlabel('Frequency (Hz)','FontSize',10);
    ylabel('Coherence','FontSize',10);
    title(sprintf('%s — peak electrode: %s', MUSCLE_LABELS{m}, peak_label),'FontSize',11);
    legend('Location','northeast','FontSize',9);
    set(ax,'FontSize',9); xlim([FREQ_MIN FREQ_MAX]); ylim([0 1]);
end

sgtitle('CMC Spectrum at Peak Beta Electrode (grand avg) — EEG vs EMG','FontSize',13);
set(gcf,'Units','pixels','Position',[50 50 1400 420]);
print(fullfile(RESULTS_DIR,'CMC_EEG_spectrum_peakElec'),'-dpng','-r300');
savefig(fullfile(RESULTS_DIR,'CMC_EEG_spectrum_peakElec.fig'));
close(gcf);

%% ===== SAVE =====
CMC_EEG.subjects          = SUBJECTS;
CMC_EEG.muscles           = MUSCLE_LABELS;
CMC_EEG.common_chanlocs   = common_chanlocs;
CMC_EEG.fr_move           = fr_move;
CMC_EEG.fr_pull           = fr_pull;
CMC_EEG.fr_base           = fr_base;
CMC_EEG.ga_move           = ga_move;    % [n_ch × n_freq_move × n_muscles]
CMC_EEG.ga_pull           = ga_pull;
CMC_EEG.ga_base           = ga_base;
CMC_EEG.band_coh_move     = band_coh_move;  % [n_ch × n_bands × n_muscles]
CMC_EEG.band_coh_pull     = band_coh_pull;
CMC_EEG.band_coh_diff     = band_coh_diff;
CMC_EEG.per_subj_move     = coh_move_subj;  % {n_subj × n_muscles}, each [n_ch × n_freq]
CMC_EEG.per_subj_pull     = coh_pull_subj;
CMC_EEG.per_subj_base     = coh_base_subj;
CMC_EEG.ci95_subj         = ci95_subj;
CMC_EEG.ci_move           = ci_move;
CMC_EEG.ci_pull           = ci_pull;
CMC_EEG.params.win_move     = WIN_MOVE;
CMC_EEG.params.win_pull     = WIN_PULL;
CMC_EEG.params.win_baseline = WIN_BASELINE;
CMC_EEG.params.bp_band      = [BP_LOW, BP_HIGH];
CMC_EEG.params.alpha        = ALPHA;

save(fullfile(RESULTS_DIR,'CMC_EEG_allSubjects.mat'), 'CMC_EEG', '-v7.3');
fprintf('\nSaved CMC_EEG results to %s\n', RESULTS_DIR);

%% ===== LOCAL FUNCTIONS =====

function [coh, fr] = mc_multitrial_coh(eeg_3d, emg_mat, fs)
% Multi-channel multi-trial Halliday coherence.
% eeg_3d:  [n_ch  × n_samp × n_trials]
% emg_mat: [n_samp × n_trials]
% coh:     [n_ch  × n_freq]  magnitude-squared coherence
% fr:      [n_freq × 1]

[n_ch, n_samp, n_trials] = size(eeg_3d);
nfft = n_samp;
win  = hann(n_samp)';     % [1 × n_samp]
nf   = nfft/2 + 1;

Sxx = zeros(n_ch, nf);   % EEG auto-spectra
Syy = zeros(1,   nf);    % EMG auto-spectrum (scalar per freq)
Sxy = zeros(n_ch, nf);   % complex cross-spectra

for t = 1:n_trials
    % EEG: [n_ch × n_samp] .* win → FFT along dim 2
    X = fft(eeg_3d(:,:,t) .* win, nfft, 2);
    X = X(:, 1:nf);                          % [n_ch × nf]

    y = emg_mat(:,t)' .* win;               % [1 × n_samp]
    Y = fft(y, nfft);
    Y = Y(1:nf);                             % [1 × nf]

    Sxx = Sxx + real(X .* conj(X));
    Syy = Syy + real(Y .* conj(Y));
    Sxy = Sxy + X .* conj(Y);               % [n_ch × nf]
end

% Halliday coherence (avoid divide-by-zero)
denom = Sxx .* Syy;
denom(denom == 0) = eps;
coh = abs(Sxy).^2 ./ denom;                 % [n_ch × nf]

fr = (0 : nf-1)' * (fs / nfft);
end


function bc = band_average(ga, fr_vec, bands)
% ga:     [n_ch × n_freq × n_muscles]
% returns [n_ch × n_bands × n_muscles]
n_ch = size(ga,1); n_mu = size(ga,3); nb = numel(bands);
bc   = zeros(n_ch, nb, n_mu);
for b = 1:nb
    fmask = fr_vec >= bands(b).lim(1) & fr_vec <= bands(b).lim(2);
    if any(fmask)
        bc(:,b,:) = mean(ga(:, fmask, :), 2);
    end
end
end


function make_topo_figure(band_coh, chanlocs, bands, muscle_labels, ...
    cond_title, cmap_name, clim_val, save_path, n_valid)
% band_coh: [n_ch × n_bands × n_muscles]
% Rows = muscles, Columns = frequency bands.

n_muscles = numel(muscle_labels);
n_bands   = numel(bands);
n_subj_v  = n_valid;

% Build colormap
switch lower(cmap_name)
    case 'hot',  cmap = hot(256);
    case 'rdbu'
        r = linspace(0,1,128)'; b_r = flip(r);
        cmap = [[b_r, b_r, ones(128,1)]; [ones(128,1), r, r]];
    otherwise,   cmap = parula(256);
end

fig = figure('Name', sprintf('CMC EEG Topoplots — %s', cond_title), ...
    'Position', [50 50 250*n_bands 250*n_muscles], 'Visible', 'off');

for m = 1:n_muscles
    for b = 1:n_bands
        ax = subplot(n_muscles, n_bands, (m-1)*n_bands + b);
        vals = band_coh(:, b, m);   % [n_ch × 1]

        topoplot(vals, chanlocs, ...
            'maplimits', clim_val, ...
            'electrodes', 'off', ...
            'colormap',   cmap, ...
            'style',      'map', ...
            'headrad',    0.5);

        if m == 1
            title(ax, bands(b).name, 'FontSize', 9);
        end
        if b == 1
            ylabel(ax, sprintf('%s\n(n=%d)', muscle_labels{m}, n_subj_v(m)), 'FontSize', 9);
        end
    end
end

sgtitle(sprintf('EEG–EMG Coherence — %s', cond_title), 'FontSize', 13);

% Single shared colorbar on the right
cb = colorbar('Position', [0.93 0.1 0.015 0.8]);
colormap(fig, cmap);
clim([clim_val(1) clim_val(2)]);
ylabel(cb, 'Coherence', 'FontSize', 9);

set(fig,'Units','pixels','Position',[50 50 250*n_bands 250*n_muscles]);
print(save_path, '-dpng', '-r200');
savefig([save_path '.fig']);
close(fig);
end
