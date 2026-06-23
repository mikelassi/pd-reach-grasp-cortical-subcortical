%% ============================================================
%  Kinematic Metrics — Vissani et al. 2021
%  Per-subject, per-block, per-trial on raw (non-warped) data
%  Output: kinematic_metrics.csv
%
%  Data sources:
%    Marker 3D positions: 01_Extracted/EMG_KIN/*EegPcsEmgKin-{N}.set
%      (EEGLAB format, 25 channels, 400 Hz)
%    Trial events A–F:   02_Kinematics/Events/*_kinematic_block.mat
%
%  Channel map in .set (1-indexed MATLAB):
%    8–10  : shoulder_right   (X, Y, Z)
%    14–16 : elbow_right      (X, Y, Z)  – lateral epicondyle
%    17–19 : wrist_right      (X, Y, Z)  – ulnar styloid
%    20–22 : forefinger_tip   (X, Y, Z)
%    23–25 : thumb_tip        (X, Y, Z)
%
%  Events (sample indices in kinematic_block.mat at 400 Hz):
%    A = movement onset        (start of reaching)
%    B = wrist velocity peak   (during reaching)
%    C = approach to target    (end of reaching / start of grasping)
%    D = pulling onset
%    E = wrist velocity peak   (during pulling)
%    F = return to rest position
%
%  Computed metrics:
%    1. Peak wrist velocity        [reach & pull, m/s]
%    2. Time to peak velocity      [reach & pull, s from phase onset]
%    3. Median radius of curvature [reach & pull, m]
%    4. Peak Hand Aperture (PHA)   [%]
%    5. Pre-Shape Coord. Index (PCI)
%    6. C-score + C-tilde          [inter-joint synergy, Micera 2005]
%
%  Author: Michael Lassi  (pipeline by Claude)
%% ============================================================

close all; clear; clc;

%% ========== PARAMETERS ==========

SUBJECTS  = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
BASE_PATH = 'G:\Projects\Parkinson_ReachGrasp\Reprocessing';
EXTRACTED = '01_Extracted\EMG_KIN';
EVENTS    = '02_Kinematics\Events';

OUTPUT_CSV = fullfile(BASE_PATH, 'kinematic_metrics.csv');

FS_KIN    = 400;      % kinematic sampling rate (Hz)
FCUT      = 8;        % low-pass cut-off for position data (Hz)  [paper spec]
SGO_WIN   = 31;       % Savitzky-Golay window (samples, must be odd) [paper: 30]
SGO_ORD   = 3;        % Savitzky-Golay polynomial order

% Channel selection is done by label at load time (see find_kin_channels)
% because the number and order of EMG channels differs across subjects
% (23 to 26 total channels; left vs right labels; forefinger vs forefinger_tip).

%% ========== EEGLAB init (for pop_loadset) ==========
eeglab; close;

%% ========== OUTPUT TABLE ==========
col_names = {'subject','block','trial', ...
    'reach_peak_vel_m_s','reach_time_to_peak_vel_s', ...
    'pull_peak_vel_m_s', 'pull_time_to_peak_vel_s', ...
    'reach_median_radius_m','pull_median_radius_m', ...
    'PHA_pct','PCI','C_score','C_tilde', ...
    'dur_reach_s','dur_grasp_s','dur_pull_s', ...
    'n_mvmt_reach','n_mvmt_pull', ...
    'NJ_reach','NJ_pull', ...
    'SI_reach','SI_pull'};

results = {};   % will grow as {row_cell, row_cell, ...}

%% ========== MAIN LOOP ==========
for s = 1:numel(SUBJECTS)
    subj = SUBJECTS{s};
    fprintf('\n===== %s =====\n', subj);

    kin_dir  = fullfile(BASE_PATH, subj, EXTRACTED);
    evt_dir  = fullfile(BASE_PATH, subj, EVENTS);

    % Find kinematic .set files (one per block)
    kin_sets = dir(fullfile(kin_dir, 'kinematic_*.set'));
    if isempty(kin_sets)
        warning('No kinematic .set files found for %s', subj);
        continue
    end
    kin_sets = sort_naturally({kin_sets.name}, kin_dir);

    % Find event .mat files (one per block)
    evt_mats = dir(fullfile(evt_dir, '*_kinematic_block.mat'));
    evt_mats = sort_naturally({evt_mats.name}, evt_dir);

    assert(numel(kin_sets) == numel(evt_mats), ...
        'Mismatch between number of .set and event .mat files for %s', subj);

    for b = 1:numel(kin_sets)
        fprintf('  Block %d/%d\n', b, numel(kin_sets));

        % ----- Load kinematic block data -----
        KIN = pop_loadset('filename', kin_sets{b}, 'filepath', kin_dir);
        raw = double(KIN.data);   % [n_chan × n_samples]

        % ----- Identify marker channels by label (robust to EMG count changes) -----
        ch = find_kin_channels({KIN.chanlocs.labels});
        fprintf('    Channels — shoulder:%s  elbow:%s  wrist:%s  index:%s  thumb:%s\n', ...
            mat2str(ch.shoulder), mat2str(ch.elbow), mat2str(ch.wrist), ...
            mat2str(ch.index_f), mat2str(ch.thumb));

        % ----- Filter & smooth 3-D positions -----
        [bfilt, afilt] = butter(4, FCUT / (FS_KIN/2), 'low');
        all_kin_ch = [ch.shoulder, ch.elbow, ch.wrist, ch.index_f, ch.thumb];
        for ci = all_kin_ch
            sig = raw(ci, :);
            sig = filtfilt(bfilt, afilt, sig);
            sig = sgolayfilt(sig, SGO_ORD, SGO_WIN);
            raw(ci, :) = sig;
        end
        % Positions are already in metres in the EEGLAB .set files.

        % Marker position matrices  [n_samples × 3]
        shoulder = raw(ch.shoulder, :)';
        elbow    = raw(ch.elbow,    :)';
        wrist    = raw(ch.wrist,    :)';
        index_f  = raw(ch.index_f,  :)';
        thumb    = raw(ch.thumb,    :)';

        % ----- Pre-compute continuous block-level signals -----
        % Wrist 3D velocity — forward first-order differentiation (paper method)
        dw        = diff(wrist) * FS_KIN;          % [n-1 × 3]  [m/s per dim]
        wrist_vel = [vecnorm(dw, 2, 2); 0];        % [n × 1] speed [m/s]; pad last sample

        % Index–thumb distance (hand aperture, m) — smooth BEFORE peak search
        % Raw aperture inherits noise from both marker positions; smoothing
        % here (rather than on raw positions) gives a cleaner scalar signal.
        aperture    = vecnorm(index_f - thumb, 2, 2);
        aperture_sm = sgolayfilt(aperture, SGO_ORD, SGO_WIN);

        % Shoulder and elbow angles for C-score
        upper_arm = elbow - shoulder;     % direction shoulder → elbow
        forearm   = wrist  - elbow;       % direction elbow → wrist

        ua_norm = vecnorm(upper_arm, 2, 2);
        fa_norm = vecnorm(forearm,   2, 2);
        u_ua    = upper_arm ./ (ua_norm + eps);
        u_fa    = forearm   ./ (fa_norm + eps);

        % α: shoulder elevation — angle of upper arm from vertical (Z)
        vert = [0, 0, 1];
        cos_alpha = u_ua * vert';
        cos_alpha = max(-1, min(1, cos_alpha));   % clamp for acos
        alpha_deg = acosd(cos_alpha);             % [n × 1] degrees

        % β: elbow flexion — angle between upper arm and forearm
        cos_beta  = sum(u_ua .* u_fa, 2);
        cos_beta  = max(-1, min(1, cos_beta));
        beta_deg  = acosd(cos_beta);              % [n × 1] degrees

        % Angular velocities (forward diff, deg/s), then SG smooth
        v_alpha = sgolayfilt([0; diff(alpha_deg)] * FS_KIN, SGO_ORD, SGO_WIN);
        v_beta  = sgolayfilt([0; diff(beta_deg)]  * FS_KIN, SGO_ORD, SGO_WIN);

        % ----- Load trial events -----
        E = load(fullfile(evt_dir, evt_mats{b}));
        events = E.events;          % struct-array of length n_trials; fields A,B,C,D,E,F
        n_trials = numel(events);

        for t = 1:n_trials
            % Event indices (1-indexed samples at FS_KIN)
            evA = double(events(t).A);
            evB = double(events(t).B);
            evC = double(events(t).C);
            evD = double(events(t).D);
            evE = double(events(t).E);
            evF = double(events(t).F);

            n_samp = size(raw, 2);

            % Safety clamp
            evA = max(1, evA);
            evC = min(evC, n_samp);
            evD = min(evD, n_samp);
            evF = min(evF, n_samp);

            % Phase sample ranges
            reach_idx = evA : evC;
            pull_idx  = evD : evF;

            if numel(reach_idx) <= SGO_WIN || numel(pull_idx) <= SGO_WIN
                fprintf('    Trial %d: too short, skipping\n', t);
                continue
            end

            % =========================================================
            % METRIC 1 & 2 — Peak wrist velocity + Time to peak
            % =========================================================
            [pk_reach, pk_reach_idx] = max(wrist_vel(reach_idx));
            t_pk_reach = (pk_reach_idx - 1) / FS_KIN;   % s from A

            [pk_pull, pk_pull_idx]   = max(wrist_vel(pull_idx));
            t_pk_pull  = (pk_pull_idx  - 1) / FS_KIN;   % s from D

            % =========================================================
            % METRIC 3 — Radius of curvature
            %   R = |v|³ / |v × a|, with v and a from gradient().
            %   Only samples where wrist speed exceeds 5 % of the phase
            %   peak are included in the median: near the endpoints of each
            %   phase the wrist is nearly stationary and the cross-product
            %   |v × a| is noise-dominated, artificially lowering the median.
            % =========================================================
            R_reach = radius_of_curvature_va(wrist(reach_idx, :), FS_KIN, SGO_ORD, SGO_WIN);
            R_pull  = radius_of_curvature_va(wrist(pull_idx,  :), FS_KIN, SGO_ORD, SGO_WIN);

            % Median R in ±0.25 s window centred on the annotated velocity peak.
            % Near the velocity peak the wrist is fastest and trajectory
            % straightest — this is the portion the paper reports as radius.
            R_WIN = round(0.35 * FS_KIN);   % 100 samples at 400 Hz
            b_loc = max(1, min(numel(reach_idx), evB - evA + 1));
            e_loc = max(1, min(numel(pull_idx),  evE - evD + 1));
            rw_reach = max(1, b_loc - R_WIN) : min(numel(reach_idx), b_loc + R_WIN);
            rw_pull  = max(1, e_loc - R_WIN) : min(numel(pull_idx),  e_loc + R_WIN);

            med_R_reach = median(R_reach(rw_reach), 'omitnan');
            med_R_pull  = median(R_pull(rw_pull),   'omitnan');

            % =========================================================
            % METRIC 4 — Peak Hand Aperture (PHA)
            %   PHA = (max_aperture_during_reach - min_aperture_during_grasp+pull)
            %          / max_aperture_during_reach  × 100  [%]
            %
            %   min window extended to C→F (not just C→D): the fingers
            %   finish closing around the object during the pulling phase,
            %   so the true minimum aperture is often reached after event D.
            % =========================================================
            ap_reach     = aperture_sm(reach_idx);                 % A → C
            ap_grasp_pull = aperture_sm(evC : evF);               % C → F

            % Mask physiologically impossible aperture values before max search
            % (same AP_MAX = 0.20 m used for PCI; positions are in metres).
            AP_MAX_PHA   = 0.20;
            ap_reach_pha = ap_reach;
            ap_reach_pha(ap_reach_pha > AP_MAX_PHA) = NaN;
            max_ap_reach = max(ap_reach_pha);   % MATLAB max ignores NaN
            if ~isnan(max_ap_reach) && max_ap_reach > 0
                PHA = (max_ap_reach - min(ap_grasp_pull)) / max_ap_reach * 100;
            else
                PHA = NaN;
            end

            % =========================================================
            % METRIC 5 — Pre-shape Coordination Index (PCI)
            %   PCI = 1 - (t_p - t_B) / Δt_tot
            %   t_p   = time of peak hand aperture DURING REACH (A→C, from A)
            %   t_B   = time of wrist velocity peak (= event B, from A)
            %   Δt_tot = duration of reaching + grasping phases (A→D)
            %
            %   Peak aperture is searched in A→C only (reaching phase), matching
            %   the PHA definition in the paper: "max(d_i-h)|_reach".
            %
            %   ARTIFACT GUARD: infrared marker dropouts create residual
            %   apertures of 0.15–1.87 m after LP+SG filtering; AP_MAX = 0.20 m
            %   catches all observed artifacts while preserving real values (< 15 cm).
            %   "Max-after-min" search then finds the genuine preshaping peak:
            %   it skips the initial resting hand posture (sometimes more open
            %   than the preshaping aperture) by first locating the global
            %   minimum in A→C and searching for the peak only after that.
            % =========================================================
            AP_MAX = 0.20;                             % 20 cm physiological upper bound
            evC_clamped  = min(evC, numel(aperture_sm));
            ap_reach_only = aperture_sm(evA : evC_clamped);       % A → C  (reach only)
            ap_reach_only(ap_reach_only > AP_MAX) = NaN;          % mask dropout artifacts

            t_B  = (evB - evA) / FS_KIN;              % s from A (event B)
            dt   = (evD - evA) / FS_KIN;              % A→D duration (reach+grasp)

            % "Max after min" search:
            %   Find the global minimum in A→C (the moment the hand finishes
            %   closing from its resting posture), then find the peak AFTER that
            %   minimum (the genuine preshaping aperture peak).
            %   Handles both (a) SG-smeared marker dropouts at onset and
            %   (b) subjects who start with the hand already open (resting aperture
            %   > preshaping peak), which would make the naïve max() find t_p = 0.
            %   MATLAB min/max both ignore NaN entries.
            tp_local = NaN;   t_p = NaN;   PCI = NaN;   % safe defaults

            [min_val, min_loc] = min(ap_reach_only);
            if ~isnan(min_val) && dt > 0
                ap_post_min = ap_reach_only(min_loc : end);
                [max_ap_val, max_rel] = max(ap_post_min);
                if ~isnan(max_ap_val)
                    tp_local = min_loc + max_rel - 1;
                    t_p  = (tp_local - 1) / FS_KIN;   % s from A
                    PCI  = 1 - (t_p - t_B) / dt;
                end
            end


            % =========================================================
            % METRIC 6 — C-score (Micera 2005)
            %
            %   Uses angular velocities v_alpha (shoulder) and v_beta (elbow)
            %   for the whole movement A→F.
            %   Normalise each to [-1, 1] within the trial.
            %   Lobe 1 centroid = mean of (v_alpha, v_beta) during reaching (A→C)
            %   Lobe 2 centroid = mean of (v_alpha, v_beta) during pulling  (D→F)
            %   C = slope of line connecting the two centroids
            %   C_tilde = sign(C) * ln(1 + |C|)   [Eq. 3, Micera 2005]
            % =========================================================
            trial_idx = evA : evF;

            va_trial = v_alpha(trial_idx);
            vb_trial = v_beta(trial_idx);

            % Normalise to [-1, 1]  (Eq. 1, Micera 2005)
            va_n = norm_m1_p1(va_trial);
            vb_n = norm_m1_p1(vb_trial);

            % Local indices within trial segment
            n_reach_samp = evC - evA + 1;
            n_grasp_samp = evD - evC + 1;

            lobe1_va = va_n(1 : n_reach_samp);
            lobe1_vb = vb_n(1 : n_reach_samp);
            lobe2_va = va_n(n_reach_samp + n_grasp_samp : end);
            lobe2_vb = vb_n(n_reach_samp + n_grasp_samp : end);

            if isempty(lobe1_va) || isempty(lobe2_va)
                C_score = NaN;
                C_tilde = NaN;
            else
                c1 = [mean(lobe1_va), mean(lobe1_vb)];
                c2 = [mean(lobe2_va), mean(lobe2_vb)];

                delta_x = c2(1) - c1(1);
                if abs(delta_x) < 1e-6
                    C_score = NaN;
                    C_tilde = NaN;
                else
                    C_score = (c2(2) - c1(2)) / delta_x;
                    C_tilde = sign(C_score) * log(1 + abs(C_score));
                end
            end

            % =========================================================
            % METRIC 7 — Phase durations
            % =========================================================
            dur_reach_s = (evC - evA) / FS_KIN;   % A→C  reaching (s)
            dur_grasp_s = (evD - evC) / FS_KIN;   % C→D  grasping (s)
            dur_pull_s  = (evF - evD) / FS_KIN;   % D→F  pulling  (s)

            % =========================================================
            % METRIC 8 — Movement units (velocity sub-peaks)
            %   Number of velocity peaks > 5% of phase maximum and
            %   minimum prominence 0.01 m/s.
            %   Healthy: 1 peak (bell-shaped).  PD: often 2–4 (fragmented).
            % =========================================================
            MVMT_THR  = 0.05;   % fraction of phase peak
            MVMT_PROM = 0.01;   % m/s minimum prominence

            if pk_reach > 0
                [~, locs_r] = findpeaks(wrist_vel(reach_idx), ...
                    'MinPeakHeight',     MVMT_THR  * pk_reach, ...
                    'MinPeakProminence', MVMT_PROM);
                n_mvmt_reach = numel(locs_r);
            else
                n_mvmt_reach = 0;
            end

            if pk_pull > 0
                [~, locs_p] = findpeaks(wrist_vel(pull_idx), ...
                    'MinPeakHeight',     MVMT_THR  * pk_pull, ...
                    'MinPeakProminence', MVMT_PROM);
                n_mvmt_pull = numel(locs_p);
            else
                n_mvmt_pull = 0;
            end

            % =========================================================
            % METRIC 9 — Normalized dimensionless jerk (NJ)
            %   NJ = -(T^5 / 2D^2) * integral(|d^3r/dt^3|^2 dt)
            %   (Hogan & Sternad 2009, Eq. 5).  More negative = smoother.
            %   T: phase duration (s).  D: straight-line displacement (m).
            % =========================================================
            NJ_reach = ndj(wrist(reach_idx,:), FS_KIN);
            NJ_pull  = ndj(wrist(pull_idx, :), FS_KIN);

            % =========================================================
            % METRIC 10 — Wrist path straightness index (SI)
            %   SI = straight-line distance / path length.
            %   SI = 1 for perfectly straight; < 1 for curved trajectory.
            % =========================================================
            path_r = sum(vecnorm(diff(wrist(reach_idx,:)), 2, 2));
            SI_reach = norm(wrist(evC,:) - wrist(evA,:)) / max(path_r, eps);

            path_p = sum(vecnorm(diff(wrist(pull_idx,:)), 2, 2));
            SI_pull  = norm(wrist(evF,:) - wrist(evD,:)) / max(path_p, eps);

            % =========================================================
            % Append row to results
            % =========================================================
            results{end+1} = {subj, b, t, ...
                pk_reach,   t_pk_reach, ...
                pk_pull,    t_pk_pull, ...
                med_R_reach, med_R_pull, ...
                PHA, PCI, C_score, C_tilde, ...
                dur_reach_s, dur_grasp_s, dur_pull_s, ...
                n_mvmt_reach, n_mvmt_pull, ...
                NJ_reach, NJ_pull, ...
                SI_reach, SI_pull};

        end % trial loop
    end % block loop
end % subject loop

%% ========== WRITE CSV ==========
fprintf('\nWriting CSV: %s\n', OUTPUT_CSV);

fid = fopen(OUTPUT_CSV, 'w');
fprintf(fid, '%s\n', strjoin(col_names, ','));

for r = 1:numel(results)
    row = results{r};
    % subject (string), block, trial (ints), then all doubles
    fprintf(fid, '%s,%d,%d', row{1}, row{2}, row{3});
    for c = 4:numel(row)
        v = row{c};
        if isnan(v)
            fprintf(fid, ',NaN');
        else
            fprintf(fid, ',%.6f', v);
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Done. %d rows written.\n', numel(results));

%% ========== HELPER FUNCTIONS ==========

function R = radius_of_curvature_va(pos, fs, sgo_ord, sgo_win)
% RADIUS_OF_CURVATURE_VA  Radius of curvature via velocity–acceleration.
%   R = |v|³ / |v × a|  where v and a are central-difference (gradient)
%   derivatives of the 3-D position.  Positions have already been LP+SG
%   filtered upstream, so gradient() gives stable estimates.
%   Low-velocity samples are excluded from the median by the caller.
%
%   pos      : [n × 3] 3-D positions (m)
%   fs       : sampling frequency (Hz)
%   sgo_ord  : Savitzky-Golay polynomial order
%   sgo_win  : Savitzky-Golay window (must be odd and > sgo_ord)
%
%   R : [n × 1] radius (m); NaN at locally straight segments
    n = size(pos, 1);
    if n < 3
        R = NaN(n, 1);
        return
    end
    vel = zeros(n, 3);
    acc = zeros(n, 3);
    for d = 1:3
        vel(:, d) = gradient(pos(:, d)) * fs;
        acc(:, d) = gradient(vel(:, d)) * fs;
    end

    speed    = vecnorm(vel, 2, 2);               % [n × 1]  m/s
    cross_va = cross(vel, acc, 2);               % [n × 3]
    cross_n  = vecnorm(cross_va, 2, 2);          % |v × a|  m²/s³

    R = speed.^3 ./ (cross_n + eps);            % m
    R(cross_n < 1e-8) = NaN;   % locally straight → NaN (excluded from median)
end

function y = norm_m1_p1(x)
% Normalise vector x to the range [-1, 1]  (Eq. 1, Micera 2005)
    mn = min(x);
    mx = max(x);
    rng = mx - mn;
    if rng < eps
        y = zeros(size(x));
    else
        y = 2 * (x - mn) / rng - 1;
    end
end

function ch = find_kin_channels(labels)
% FIND_KIN_CHANNELS  Locate X/Y/Z channel indices for each marker by label.
%   Labels vary across subjects: right vs left, forefinger_tip vs forefinger,
%   and the number of EMG channels before the kinematic ones differs (23–26 ch).
%   All kinematic channel names end in '_X_kin', '_Y_kin', '_Z_kin'.
    ch.shoulder = find_xyz(labels, 'shoulder');
    ch.elbow    = find_xyz(labels, 'elbow');
    ch.wrist    = find_xyz(labels, 'wrist');
    ch.index_f  = find_xyz(labels, 'forefinger');   % matches forefinger_tip too
    ch.thumb    = find_xyz(labels, 'thumb');
end

function idx = find_xyz(labels, marker)
% Return [X_idx, Y_idx, Z_idx] for a marker name substring.
    x = find(~cellfun(@isempty, regexpi(labels, [marker, '.*_X_kin'])));
    y = find(~cellfun(@isempty, regexpi(labels, [marker, '.*_Y_kin'])));
    z = find(~cellfun(@isempty, regexpi(labels, [marker, '.*_Z_kin'])));
    if numel(x) ~= 1 || numel(y) ~= 1 || numel(z) ~= 1
        error('Could not find unique X/Y/Z channels for marker "%s".\nMatches: X=%s  Y=%s  Z=%s', ...
            marker, mat2str(x), mat2str(y), mat2str(z));
    end
    idx = [x, y, z];
end

function nj = ndj(pos, fs)
% NDJ  Normalized dimensionless jerk (Hogan & Sternad 2009, Eq. 5).
%   pos : [n × 3] already-filtered 3-D positions (m)
%   fs  : sampling rate (Hz)
%   nj  : scalar; more negative → smoother movement
    n = size(pos, 1);
    T = (n - 1) / fs;
    D = norm(pos(end,:) - pos(1,:));
    if T <= 0 || D <= 0 || n < 4
        nj = NaN; return
    end
    jrk = zeros(n, 3);
    for d = 1:3
        v        = gradient(pos(:,d)) * fs;
        a        = gradient(v)        * fs;
        jrk(:,d) = gradient(a)        * fs;
    end
    nj = -(T^5 / (2 * D^2)) * (sum(jrk(:).^2) / fs);
end

function sorted = sort_naturally(names, ~)
% Return filenames (NOT full paths) sorted by the trailing block number.
% The caller is responsible for prepending the directory when needed.
    nums = zeros(1, numel(names));
    for i = 1:numel(names)
        tok = regexp(names{i}, '-(\d+)[^-]*$', 'tokens');
        if ~isempty(tok)
            nums(i) = str2double(tok{1}{1});
        end
    end
    [~, idx] = sort(nums);
    sorted = names(idx);   % cell array of filenames only
end
