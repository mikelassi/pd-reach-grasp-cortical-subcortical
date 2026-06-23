%% PreEpoching.m
% =========================================================================
% Pre-epoching analysis for reach-to-grasp LFP data
%
% Purpose:
%   Load unsegmented, preprocessed LFP files (with kinematic events A-F),
%   compute the duration of each kinematic phase across all trials, blocks,
%   and subjects, then recommend conservative epoch lengths for downstream
%   analysis.
%
% Kinematic event definitions (from A02_events_segmentation.m):
%   A : movement onset     (wrist velocity crosses threshold → reach start)
%   B : peak reach velocity
%   C : velocity minimum   (reach end / object contact → grasp start)
%   D : pull onset         (velocity re-crosses threshold → grasp end)
%   E : peak pull velocity
%   F : return to baseline (end of pull)
%
% Phases used here:
%   Reach  : A → C
%   Grasp  : C → D
%   Pull   : D → F
%   ITI    : F_t → A_{t+1}  (inter-trial interval = available rest)
%             [only interior intervals; edge trials excluded for rest]
%
% Output:
%   - Console table: min / p5 / p10 / p25 / p50 / p75 / max per phase
%   - Box plots of phase duration distributions with recommended cutoff
%   - Empirical CDFs per phase
%   - Schematic of a single trial showing phase boundaries
%   - Printed recommendation block (copy-paste into epoching script)
%
% Author: Michael Lassi
% =========================================================================

close all; clear; clc;

%% ===== PARAMETERS =====

SUBJECTS        = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH       = 'H:\Parkinson_ReachGrasp\Reprocessing';
PREPROC_SUB     = fullfile('Preprocessed', 'LFP');

% Percentile for conservative epoch-length recommendation.
%   5  → 95 % of all trials fit entirely within the epoch
%   10 → 90 % of all trials fit
EPOCH_PERCENTILE = 5;

% Minimum gap (s) to include as rest padding around movement epochs.
% Used only for the recommendation summary; does not filter trials.
MIN_REST_BUFFER_S = 0.5;

%% ===== INIT =====
eeglab; close;

% Pre-allocate storage (will grow; ~300 trials expected)
% Columns: [reach, grasp, pull, iti]  (samples)
all_dur   = [];   % N_total × 4
meta      = {};   % N_total × 3  {subject, block, trial}

fs = NaN;

PHASE_IDX  = struct('reach', 1, 'grasp', 2, 'pull', 3, 'iti', 4);
PHASE_COLS = {[0.3 0.8 0.3], [0.2 0.6 1.0], [0.6 0.1 0.9], [0.8 0.2 0.2]};
PHASE_NAMES = {'Reach (A→C)', 'Grasp (C→D)', 'Pull (D→F)', 'ITI (F→A_{next})'};

%% ===== MAIN LOOP =====
for s = 1:numel(SUBJECTS)
    subj         = SUBJECTS{s};
    preproc_path = fullfile(BASE_PATH, subj, PREPROC_SUB);
    set_files    = dir(fullfile(preproc_path, '*_wEv.set'));

    if isempty(set_files)
        warning('No *_wEv.set files found for %s in:\n  %s', subj, preproc_path);
        continue;
    end

    fprintf('\n===== Subject: %s  (%d block(s)) =====\n', subj, numel(set_files));

    for f = 1:numel(set_files)
        EEG = pop_loadset('filename', set_files(f).name, 'filepath', preproc_path);

        if isnan(fs)
            fs = EEG.srate;
            fprintf('  Sampling rate detected: %d Hz\n', fs);
        end

        % ------------------------------------------------------------------
        % Parse events: keep only those matching the pattern [A-F]_T<number>
        % (ignores BrainVision / TENStrigger / other markers in the file)
        % ------------------------------------------------------------------
        ev_types    = {EEG.event.type};
        ev_lats     = [EEG.event.latency];   % in samples

        trial_nums = [];
        for e = 1:numel(ev_types)
            tok = regexp(ev_types{e}, '^[A-F]_T(\d+)$', 'tokens');
            if ~isempty(tok)
                trial_nums(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
            end
        end

        if isempty(trial_nums)
            warning('  No kinematic events in block %d of %s — skipping.', f, subj);
            continue;
        end

        trial_nums = unique(trial_nums);
        n_trials   = numel(trial_nums);

        % Build latency matrix: rows = trials, cols = [A B C D E F]
        EV_LETTERS = 'ABCDEF';
        lat = nan(n_trials, 6);
        for e = 1:numel(ev_types)
            tok = regexp(ev_types{e}, '^([A-F])_T(\d+)$', 'tokens');
            if isempty(tok); continue; end
            letter = tok{1}{1};
            t_num  = str2double(tok{1}{2});
            t_idx  = find(trial_nums == t_num);
            e_idx  = find(EV_LETTERS == letter);
            if ~isempty(t_idx) && ~isempty(e_idx)
                lat(t_idx, e_idx) = ev_lats(e);
            end
        end

        % ------------------------------------------------------------------
        % Compute phase durations for each trial
        % ITI is only available for interior boundaries (not first/last)
        % ------------------------------------------------------------------
        n_parsed = 0;
        for t = 1:n_trials
            A = lat(t, 1);  C = lat(t, 3);  D = lat(t, 4);  F = lat(t, 6);

            if any(isnan([A C D F]))
                fprintf('  [SKIP] %s blk%d trial %d — missing A/C/D/F\n', ...
                    subj, f, trial_nums(t));
                continue;
            end

            reach = C - A;
            grasp = D - C;
            pull  = F - D;

            % Sanity: all durations must be positive
            if any([reach, grasp, pull] <= 0)
                fprintf('  [WARN] %s blk%d trial %d — non-positive phase (R=%d G=%d P=%d)\n', ...
                    subj, f, trial_nums(t), reach, grasp, pull);
                continue;
            end

            % ITI: only for trials where the next trial's A is available
            if t < n_trials
                A_next = lat(t+1, 1);
                iti    = A_next - F;
                if isnan(A_next) || iti <= 0
                    iti = NaN;
                end
            else
                iti = NaN;   % last trial in block — no following trial
            end

            all_dur(end+1, :) = [reach, grasp, pull, iti]; %#ok<AGROW>
            meta(end+1, :)    = {subj, f, trial_nums(t)};   %#ok<AGROW>
            n_parsed = n_parsed + 1;
        end

        fprintf('  Block %d: %d / %d trials parsed\n', f, n_parsed, n_trials);
    end
end

assert(~isnan(fs), 'No files were loaded — check BASE_PATH and PREPROC_SUB.');
n_total = size(all_dur, 1);
fprintf('\nTotal trials parsed: %d\n', n_total);

%% ===== STATISTICS TABLE =====
dur_s  = all_dur / fs;   % convert to seconds

PCT_LEVELS = [0, 5, 10, 25, 50, 75, 90, 95, 100];
n_phases   = numel(PHASE_NAMES);

fprintf('\n');
fprintf('================================================================\n');
fprintf('  PHASE DURATION STATISTICS  (n = %d trials, fs = %d Hz)\n', n_total, fs);
fprintf('================================================================\n');

hdr_fmt = '%-20s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s\n';
row_fmt = '%-20s  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f\n';
fprintf(hdr_fmt, 'Phase', 'min', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'max');
fprintf('%s\n', repmat('-', 1, 68));

pct_table = nan(n_phases, numel(PCT_LEVELS));
n_valid   = zeros(1, n_phases);

for p = 1:n_phases
    vals = dur_s(:, p);
    vals = vals(~isnan(vals));
    n_valid(p) = numel(vals);
    pct_table(p, :) = prctile(vals, PCT_LEVELS);
    fprintf(row_fmt, ...
        [PHASE_NAMES{p} sprintf(' (n=%d)', n_valid(p))], pct_table(p, :));
end

%% ===== RECOMMENDATION =====
fprintf('\n');
fprintf('================================================================\n');
fprintf('  RECOMMENDED EPOCH LENGTHS  (%dth pct → %.0f%% of trials fit)\n', ...
    EPOCH_PERCENTILE, 100 - EPOCH_PERCENTILE);
fprintf('================================================================\n');
fprintf('%-20s  %8s  %8s  %8s\n', 'Phase', 'Seconds', 'Samples', 'Covers');
fprintf('%s\n', repmat('-', 1, 52));

rec_s    = zeros(1, n_phases);
rec_samp = zeros(1, n_phases);

for p = 1:n_phases
    vals          = dur_s(:, p);
    vals          = vals(~isnan(vals));
    rec_s(p)      = prctile(vals, EPOCH_PERCENTILE);
    rec_samp(p)   = floor(rec_s(p) * fs);
    pct_covered   = mean(vals <= rec_s(p)) * 100;
    fprintf('%-20s  %8.3f  %8d  %6.1f%%\n', ...
        PHASE_NAMES{p}, rec_s(p), rec_samp(p), pct_covered);
end

% Available ITI buffer (subtract recommended ITI epoch to find usable rest margin)
% ITI is shared between rest-after of trial t and rest-before of trial t+1.
% Each side gets half of the ITI (conservative).
iti_vals      = dur_s(~isnan(dur_s(:, PHASE_IDX.iti)), PHASE_IDX.iti);
iti_p5        = prctile(iti_vals, EPOCH_PERCENTILE);
iti_half_p5   = iti_p5 / 2;

fprintf('\n');
fprintf('Conservative ITI (%dth pct) : %.3f s → each side: %.3f s\n', ...
    EPOCH_PERCENTILE, iti_p5, iti_half_p5);
fprintf('Requested MIN_REST_BUFFER   : %.3f s\n', MIN_REST_BUFFER_S);

if iti_half_p5 >= MIN_REST_BUFFER_S
    fprintf('  → Buffer FEASIBLE: %.3f s available per side.\n', iti_half_p5);
else
    fprintf('  → WARNING: Only %.3f s available per side (< %.3f s requested).\n', ...
        iti_half_p5, MIN_REST_BUFFER_S);
    fprintf('     Consider reducing MIN_REST_BUFFER_S or reducing EPOCH_PERCENTILE.\n');
end

% --- Copy-paste summary ---
fprintf('\n');
fprintf('================================================================\n');
fprintf('  COPY-PASTE PARAMETERS FOR EPOCHING SCRIPT\n');
fprintf('================================================================\n');
fprintf('EPOCH_REACH_S  = %.3f;   %% %dth-pct reach duration (s)\n',  rec_s(1), EPOCH_PERCENTILE);
fprintf('EPOCH_GRASP_S  = %.3f;   %% %dth-pct grasp duration (s)\n',  rec_s(2), EPOCH_PERCENTILE);
fprintf('EPOCH_PULL_S   = %.3f;   %% %dth-pct pull duration (s)\n',   rec_s(3), EPOCH_PERCENTILE);
fprintf('EPOCH_REST_S   = %.3f;   %% half ITI, %dth pct (s)\n', iti_half_p5, EPOCH_PERCENTILE);

%% ===== PLOT 1: Box plots of phase durations =====
figure('Name', 'Phase Duration Distributions', 'Position', [50 50 1400 500]);

for p = 1:n_phases
    subplot(1, n_phases, p);

    vals = dur_s(:, p);
    vals = vals(~isnan(vals));

    bh = boxplot(vals, 'Colors', PHASE_COLS{p}, 'Symbol', '+k', 'Widths', 0.5);
    set(bh, 'LineWidth', 1.5);
    hold on;

    % Mark recommended cutoff
    yline(rec_s(p), '--k', sprintf('p%d: %.2fs', EPOCH_PERCENTILE, rec_s(p)), ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment',   'bottom', ...
        'LineWidth', 1.8, 'FontSize', 9);

    % Mark median
    yline(prctile(vals, 50), ':k', sprintf('p50: %.2fs', prctile(vals,50)), ...
        'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment',   'top', ...
        'LineWidth', 1.2, 'FontSize', 8);

    title(PHASE_NAMES{p}, 'FontSize', 11);
    ylabel('Duration (s)', 'FontSize', 10);
    xlabel(sprintf('n = %d', n_valid(p)), 'FontSize', 9);
    set(gca, 'XTick', [], 'FontSize', 10);
    grid on;
end

sgtitle(sprintf('Phase Duration Distributions — All Subjects & Trials  (dashed = p%d)', ...
    EPOCH_PERCENTILE), 'FontSize', 13);

%% ===== PLOT 2: Empirical CDFs =====
figure('Name', 'Cumulative Phase Duration Distributions', 'Position', [50 600 900 420]);
hold on; grid on;

for p = 1:n_phases
    vals = sort(dur_s(:, p));
    vals = vals(~isnan(vals));
    cdf  = (1:numel(vals))' / numel(vals) * 100;
    plot(vals, cdf, 'Color', PHASE_COLS{p}, 'LineWidth', 2.2, 'DisplayName', PHASE_NAMES{p});
end

% Horizontal reference lines at key percentiles
for pct_ref = [50, 90, 95]
    yline(pct_ref, ':k', sprintf('%d%%', pct_ref), ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8, 'LineWidth', 0.8);
end

xlabel('Duration (s)', 'FontSize', 12);
ylabel('Trials fitting within duration (%)', 'FontSize', 12);
title('Empirical CDFs of Phase Durations', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 10);
set(gca, 'FontSize', 11);
xlim([0, prctile(dur_s(:), 99)]);

%% ===== PLOT 3: Per-subject median durations (check for outlier subjects) =====
figure('Name', 'Per-Subject Median Phase Durations', 'Position', [1000 50 700 500]);
hold on; grid on;

subj_meds = nan(numel(SUBJECTS), 4);
for s = 1:numel(SUBJECTS)
    mask = strcmp(meta(:,1), SUBJECTS{s});
    if ~any(mask); continue; end
    for p = 1:4
        vals = dur_s(mask, p);
        subj_meds(s, p) = median(vals(~isnan(vals)));
    end
end

x_pos = 1:numel(SUBJECTS);
for p = 1:4
    plot(x_pos, subj_meds(:, p), '-o', ...
        'Color', PHASE_COLS{p}, 'LineWidth', 2, 'MarkerSize', 8, ...
        'MarkerFaceColor', PHASE_COLS{p}, 'DisplayName', PHASE_NAMES{p});
end

set(gca, 'XTick', x_pos, 'XTickLabel', SUBJECTS, 'XTickLabelRotation', 30);
xlabel('Subject', 'FontSize', 12);
ylabel('Median duration (s)', 'FontSize', 12);
title('Per-Subject Median Phase Durations', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

%% ===== PLOT 4: Schematic trial timeline =====
figure('Name', 'Trial Phase Schematic', 'Position', [1000 600 900 250]);
ax = axes; hold(ax, 'on');

% Median durations for schematic
med_reach = prctile(dur_s(~isnan(dur_s(:,1)),1), 50);
med_grasp = prctile(dur_s(~isnan(dur_s(:,2)),2), 50);
med_pull  = prctile(dur_s(~isnan(dur_s(:,3)),3), 50);
med_iti   = prctile(iti_vals, 50);

t_A = med_iti / 2;
t_C = t_A + med_reach;
t_D = t_C + med_grasp;
t_F = t_D + med_pull;
t_end = t_F + med_iti / 2;

phase_starts  = [0,    t_A,  t_C,  t_D,  t_F];
phase_ends    = [t_A,  t_C,  t_D,  t_F,  t_end];
phase_labels  = {'Rest', 'Reach', 'Grasp', 'Pull', 'Rest'};
phase_colours = {[0.8 0.2 0.2], [0.3 0.8 0.3], [0.2 0.6 1.0], [0.6 0.1 0.9], [1.0 0.55 0.1]};

for p = 1:5
    rectangle(ax, 'Position', [phase_starts(p), 0, phase_ends(p)-phase_starts(p), 1], ...
        'FaceColor', [phase_colours{p}, 0.5], 'EdgeColor', 'none');
    text(ax, (phase_starts(p)+phase_ends(p))/2, 0.5, phase_labels{p}, ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end

for xv = [t_A, t_C, t_D, t_F]
    xline(ax, xv, 'k-', 'LineWidth', 2);
end

for xv_lab = {t_A, 'A'; t_C, 'C'; t_D, 'D'; t_F, 'F'}'
    text(ax, xv_lab{1}, -0.12, xv_lab{2}, 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Color', 'k');
end

xlim([0, t_end]); ylim([-0.3, 1.2]);
ax.YAxis.Visible = 'off';
ax.XTick = [t_A, t_C, t_D, t_F];
ax.XTickLabel = arrayfun(@(v) sprintf('%.2fs', v), [t_A, t_C, t_D, t_F], ...
    'UniformOutput', false);
xlabel('Time (s) — median durations', 'FontSize', 11);
title('Trial Phase Structure (Schematic, Median Durations)', 'FontSize', 12);
set(ax, 'FontSize', 10);
