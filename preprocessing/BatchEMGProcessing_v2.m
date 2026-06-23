%% BatchEMGProcessing_v2.m
% =========================================================================
% Batch EMG preprocessing for all 8 PD reach-to-grasp subjects.
%
% Input:
%   H:\Parkinson_ReachGrasp\Reprocessing\<subj>\01_Extracted\EMG_KIN\*.set
%
% Channel selection (label-based, not index-based):
%   - wue02, wue03, wue10, wue11 have both dominant + non-dominant EMG sets.
%     Contralateral to LFP (right STN / left hand movement = LEFT arm)
%     → select channels with '_dominant_emg' in their label.
%   - wue05, wue06, wue07, wue09 have a single EMG set.
%     → select all channels containing '_emg' (case-insensitive).
%
%   After contralateral selection, restrict to the 3 muscles common to all
%   8 subjects: IOD, Triceps, Deltoid.
%
% Pipeline (per block, per subject):
%   1. Demean
%   2. Band-pass  20–180 Hz  (4th-order Butterworth, zero-phase)
%   3. Notch      50/100 Hz  (2nd-order Butterworth notch, zero-phase)
%   4. Rectify    (absolute value of band-passed signal)
%   5. LP envelope 5 Hz      (4th-order Butterworth, zero-phase)
%   6. Normalise  to 95th-percentile of the envelope (per channel)
%
% Outputs (saved to RESULTS_DIR/EMG/<subj>/):
%   emg_bp_<subj>.set      — band-passed signal (step 2 only)
%   emg_rect_<subj>.set    — rectified signal   (steps 1–4)
%   emg_env_<subj>.set     — normalised envelope (all steps)
%   EMG_allSubjects.mat    — struct with envelope, rect, bp + metadata
%
% Author: Michael Lassi
% =========================================================================

clear; clc;

%% ===== PARAMETERS =====
SUBJECTS    = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH   = 'H:\Parkinson_ReachGrasp\Reprocessing';
INPUT_SUB   = fullfile('01_Extracted', 'EMG_KIN');
RESULTS_DIR = fullfile(BASE_PATH, 'RESULTS_final', 'EMG');

% Subjects with separate dominant / non-dominant EMG sets.
% Contralateral = dominant (right STN / left-arm movers).
DOMINANT_SUBJECTS = {'wue02','wue03','wue10','wue11'};

% The 3 muscles present in every subject (case-insensitive partial match).
% Script will warn if any are absent for a given subject.
TARGET_MUSCLES = {'iod', 'triceps', 'deltoid'};

% Pipeline parameters
BP_LOW   = 20;      % Hz  — band-pass lower edge
BP_HIGH  = 180;     % Hz  — band-pass upper edge
BP_ORDER = 4;
NOTCH_FREQS = [50, 100];   % Hz
NOTCH_BW    = 2;           % ± Hz half-bandwidth for notch
NOTCH_ORDER = 2;
ENV_CUTOFF  = 5;    % Hz  — low-pass cutoff for envelope
ENV_ORDER   = 4;
NORM_PCTILE = 95;   % percentile for amplitude normalisation

%% ===== BIOSIG STUB FIX =====
% EEGLAB's Biosig plugin ships stub versions of butter/filtfilt that shadow
% the Signal Processing Toolbox equivalents — remove them before filtering.
biosig_stubs = '';
if ~isempty(which('eeglab'))
    biosig_stubs = fullfile(fileparts(which('eeglab')), ...
        'plugins', 'Biosig3.8.4', 'biosig', 'maybe-missing');
end
if ~isempty(biosig_stubs) && exist(biosig_stubs, 'dir')
    rmpath(biosig_stubs);
    fprintf('Removed Biosig stubs from path.\n');
end

%% ===== INIT =====
eeglab nogui;

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

EMG_ALL = struct();   % consolidated output struct

proc_log = {};        % processing log (cell array of strings)

%% ===== SUBJECT LOOP =====
for s = 1:numel(SUBJECTS)
    subj      = SUBJECTS{s};
    emg_path  = fullfile(BASE_PATH, subj, INPUT_SUB);
    set_files = dir(fullfile(emg_path, '*.set'));

    if isempty(set_files)
        proc_log{end+1} = sprintf('[WARN] %s: no .set files in %s', subj, emg_path); %#ok<AGROW>
        warning('No .set files found for %s — skipping.', subj);
        continue;
    end

    fprintf('\n===== Subject: %s  (%d block(s)) =====\n', subj, numel(set_files));

    % Determine whether to use dominant or all EMG labels
    use_dominant = ismember(subj, DOMINANT_SUBJECTS);

    % Accumulators across blocks
    EEG_bp_all   = [];
    EEG_rect_all = [];
    EEG_env_all  = [];

    for f = 1:numel(set_files)
        EEG = pop_loadset('filename', set_files(f).name, 'filepath', emg_path);
        fs  = EEG.srate;
        fprintf('  Block %d | %s | %d ch, %d samp\n', f, set_files(f).name, EEG.nbchan, EEG.pnts);

        % ---------------------------------------------------------------
        % STEP 1: Select contralateral EMG channels by label
        % ---------------------------------------------------------------
        ch_labels = lower({EEG.chanlocs.labels});

        if use_dominant
            sel_mask = contains(ch_labels, '_dominant_emg');
        else
            sel_mask = contains(ch_labels, '_emg');
        end

        if ~any(sel_mask)
            proc_log{end+1} = sprintf('[WARN] %s block %d: no matching EMG channels found', subj, f); %#ok<AGROW>
            warning('  Block %d: no matching EMG channels for %s — skipping block.', f, subj);
            continue;
        end

        ch_idx = find(sel_mask);

        % Restrict to TARGET_MUSCLES (IOD, Triceps, Deltoid)
        target_mask = false(1, numel(ch_idx));
        for m = 1:numel(TARGET_MUSCLES)
            hits = contains(ch_labels(ch_idx), TARGET_MUSCLES{m});
            if ~any(hits)
                proc_log{end+1} = sprintf('[WARN] %s block %d: muscle "%s" not found', ...
                    subj, f, TARGET_MUSCLES{m}); %#ok<AGROW>
                warning('  Block %d: muscle "%s" not found for %s.', f, TARGET_MUSCLES{m}, subj);
            end
            target_mask = target_mask | hits;
        end
        ch_idx = ch_idx(target_mask);

        if isempty(ch_idx)
            proc_log{end+1} = sprintf('[WARN] %s block %d: none of target muscles found', subj, f); %#ok<AGROW>
            warning('  Block %d: none of target muscles found for %s — skipping.', f, subj);
            continue;
        end

        sel_labels = {EEG.chanlocs(ch_idx).labels};
        fprintf('    Selected %d channels: %s\n', numel(ch_idx), strjoin(sel_labels, ', '));

        % Keep only selected channels
        EEG_sel = pop_select(EEG, 'channel', ch_idx);

        % ---------------------------------------------------------------
        % STEP 2: Demean
        % ---------------------------------------------------------------
        raw = double(EEG_sel.data);   % [n_ch × n_samp]
        raw = raw - mean(raw, 2);

        % ---------------------------------------------------------------
        % STEP 3: Band-pass 20–180 Hz
        % ---------------------------------------------------------------
        [b_bp, a_bp] = butter(BP_ORDER, [BP_LOW, BP_HIGH] / (fs/2), 'bandpass');
        bp = zeros(size(raw));
        for c = 1:size(raw, 1)
            bp(c,:) = filtfilt(b_bp, a_bp, raw(c,:));
        end

        % ---------------------------------------------------------------
        % STEP 4: Notch at 50 Hz and harmonics
        % ---------------------------------------------------------------
        sig_notched = bp;
        for fn = NOTCH_FREQS
            if fn < fs/2
                bw = NOTCH_BW / (fs/2);
                fc = fn / (fs/2);
                [b_n, a_n] = butter(NOTCH_ORDER, [fc - bw/2, fc + bw/2], 'stop');
                for c = 1:size(sig_notched, 1)
                    sig_notched(c,:) = filtfilt(b_n, a_n, sig_notched(c,:));
                end
            end
        end

        % ---------------------------------------------------------------
        % STEP 5: Rectify
        % ---------------------------------------------------------------
        rect = abs(sig_notched);

        % ---------------------------------------------------------------
        % STEP 6: LP envelope at 5 Hz
        % ---------------------------------------------------------------
        [b_env, a_env] = butter(ENV_ORDER, ENV_CUTOFF / (fs/2), 'low');
        env = zeros(size(rect));
        for c = 1:size(rect, 1)
            env(c,:) = filtfilt(b_env, a_env, rect(c,:));
        end

        % ---------------------------------------------------------------
        % STEP 7: Normalise to 95th percentile (per channel)
        % ---------------------------------------------------------------
        env_norm = zeros(size(env));
        for c = 1:size(env, 1)
            p95 = prctile(env(c,:), NORM_PCTILE);
            if p95 > 0
                env_norm(c,:) = env(c,:) / p95;
            else
                env_norm(c,:) = env(c,:);
                proc_log{end+1} = sprintf('[WARN] %s block %d ch %d: 95th pct=0, normalisation skipped', ...
                    subj, f, c); %#ok<AGROW>
            end
        end

        % ---------------------------------------------------------------
        % Pack back into EEGLAB structs
        % ---------------------------------------------------------------
        EEG_bp_f        = EEG_sel; EEG_bp_f.data   = single(sig_notched);
        EEG_rect_f      = EEG_sel; EEG_rect_f.data = single(rect);
        EEG_env_f       = EEG_sel; EEG_env_f.data  = single(env_norm);

        % Merge across blocks
        if isempty(EEG_bp_all)
            EEG_bp_all   = EEG_bp_f;
            EEG_rect_all = EEG_rect_f;
            EEG_env_all  = EEG_env_f;
        else
            EEG_bp_all   = pop_mergeset(EEG_bp_all,   EEG_bp_f);
            EEG_rect_all = pop_mergeset(EEG_rect_all, EEG_rect_f);
            EEG_env_all  = pop_mergeset(EEG_env_all,  EEG_env_f);
        end

    end  % block loop

    if isempty(EEG_bp_all)
        proc_log{end+1} = sprintf('[SKIP] %s: no valid data after channel selection', subj); %#ok<AGROW>
        warning('Subject %s: no valid EMG data — skipping save.', subj);
        continue;
    end

    % ---------------------------------------------------------------
    % Save per-subject .set files
    % ---------------------------------------------------------------
    subj_out = fullfile(RESULTS_DIR, subj);
    if ~exist(subj_out, 'dir'), mkdir(subj_out); end

    pop_saveset(EEG_bp_all,   'filename', sprintf('emg_bp_%s.set',   subj), 'filepath', subj_out);
    pop_saveset(EEG_rect_all, 'filename', sprintf('emg_rect_%s.set', subj), 'filepath', subj_out);
    pop_saveset(EEG_env_all,  'filename', sprintf('emg_env_%s.set',  subj), 'filepath', subj_out);
    fprintf('  Saved: emg_bp/rect/env_%s.set  (%d total samples, %d ch)\n', ...
        subj, EEG_env_all.pnts, EEG_env_all.nbchan);

    % ---------------------------------------------------------------
    % Store in consolidated struct
    % ---------------------------------------------------------------
    EMG_ALL.(subj).labels   = {EEG_env_all.chanlocs.labels};
    EMG_ALL.(subj).fs       = EEG_env_all.srate;
    EMG_ALL.(subj).envelope = double(EEG_env_all.data);   % [n_ch × n_samp]
    EMG_ALL.(subj).rect     = double(EEG_rect_all.data);
    EMG_ALL.(subj).bp       = double(EEG_bp_all.data);

    % ---------------------------------------------------------------
    % Diagnostic plot: mean envelope per channel
    % ---------------------------------------------------------------
    fig = figure('Name', sprintf('EMG Envelope — %s', subj), ...
        'Position', [50 50 1200 400], 'Visible', 'off');
    n_ch   = EEG_env_all.nbchan;
    t_axis = (0 : EEG_env_all.pnts - 1) / EEG_env_all.srate;
    for c = 1:n_ch
        subplot(1, n_ch, c);
        plot(t_axis, EMG_ALL.(subj).envelope(c,:), 'LineWidth', 0.5);
        title(EMG_ALL.(subj).labels{c}, 'Interpreter', 'none', 'FontSize', 9);
        xlabel('Time (s)'); ylabel('Norm. amplitude');
        ylim([0 3]);
    end
    sgtitle(sprintf('%s — EMG Normalised Envelope', subj), 'FontSize', 12);
    print(fullfile(subj_out, sprintf('EMG_envelope_%s', subj)), '-dpng', '-r150');
    close(fig);

end  % subject loop

%% ===== SAVE CONSOLIDATED .mat =====
EMG_ALL.params.subjects     = SUBJECTS;
EMG_ALL.params.target_muscles = TARGET_MUSCLES;
EMG_ALL.params.bp_band      = [BP_LOW, BP_HIGH];
EMG_ALL.params.env_cutoff   = ENV_CUTOFF;
EMG_ALL.params.norm_pctile  = NORM_PCTILE;
EMG_ALL.params.proc_log     = proc_log;

save(fullfile(RESULTS_DIR, 'EMG_allSubjects.mat'), 'EMG_ALL', '-v7.3');
fprintf('\nSaved consolidated EMG results to %s\n', RESULTS_DIR);

%% ===== PRINT LOG =====
if ~isempty(proc_log)
    fprintf('\n--- Processing log ---\n');
    for k = 1:numel(proc_log)
        fprintf('%s\n', proc_log{k});
    end
end
fprintf('\nDone.\n');
