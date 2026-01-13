% =========================================================================
% Script for automatic detection of kinematic events on EMG/KIN data
% =========================================================================
% GENERAL DESCRIPTION:
% This code processes EEG/EMG/KIN files from a list of subjects,
% automatically extracting and classifying six kinematic events
% (A, B, C, D, E, F) from the absolute wrist velocity signal.
%
% STRUCTURE AND WORKFLOW:
% 1. Sets subject-specific parameters (search windows, thresholds, lags).
% 2. Loads the .set files for each subject from the "Extracted/EMG_KIN" folder.
% 3. Identifies wrist kinematic channels (left or right, depending on the subject).
% 4. Applies a low-pass filter (8 Hz) and Savitzky–Golay smoothing (31 samples).
% 5. Computes absolute wrist velocity (temporal derivative and Euclidean norm).
% 6. Applies additional smoothing (moving average ~100 ms) to reduce residual noise.
% 7. Determines an initial threshold based on baseline activity
%    (mean + k * standard deviation).
% 8. Scans the data to identify the temporal locations of six events:
%       - A: movement onset exceeding the threshold
%       - B: maximum velocity peak
%       - C: first minimum below a threshold after B
%       - D: re-onset exceeding a local threshold
%       - E: new velocity peak
%       - F: return to baseline velocity
% 9. Plots the velocity signal with color-coded markers for each event.
% 10. Saves high-resolution PNG figures in the subject folder.
% 11. Saves the "events" structure to a .mat file, unless manual events
%     are required (special cases → loaded and plotted instead of automatic ones).
%
% NOTES:
% - Parameters (lags, thresholds, offsets) are tuned on a subject-by-subject basis.
% - The baseline window is defined in milliseconds, while latencies are in samples.
% - Pre-saved manual events are loaded and plotted when required.
% - The code also generates plots for visual inspection of the results.
% =========================================================================

close all;

subject_list = {'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'}; 
base_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi\';
eeglab;


for s = 1:length(subject_list)
    subject_id = subject_list{s};

    fprintf('=== Processing subject %s ===\n', subject_id);
     
    % Ad hoc parameters
    % Offset for D event threshold
    if any(strcmp(subject_id, {'wue05', 'wue03'})) 
        D_add = 0.03;
    else
        D_add = 0.01;
    end
    
    % Search restart parameter after F identification
    if any(strcmp(subject_id, {'wue02'})) 
        F_lag = 100;
    elseif any(strcmp(subject_id, {'wue05'})) 
        F_lag = 400;
    else
        F_lag = 200;
    end
    
    if any(strcmp(subject_id, {'wue05'})) % Start search 1000 samples later
        start_offset = 1000;
    elseif any(strcmp(subject_id, {'wue07'})) % Start search 2600 samples later
        start_offset = 2600;
    else
        start_offset = 0;
    end


    % EMG_KIN folder path
    set_path = fullfile(base_path, subject_id, 'Extracted', 'EMG_KIN');

    % List of .set files in the folder
    set_files = dir(fullfile(set_path, '*.set'));

    for f = 1:length(set_files)
        set_file = set_files(f).name;
        fprintf('  > File: %s\n', set_file);

        % Dataset loading
        EEG = pop_loadset('filename', set_file, 'filepath', set_path);
        time = EEG.times;

        % Wrist channel extraction (for some subjects the left wrist
        % is available, for others the right one)
        if any(strcmp(subject_id, {'wue05', 'wue09'}))
            ch_wrist = find(ismember({EEG.chanlocs.labels}, ...
                {'wrist_left_X_kin', 'wrist_left_Y_kin', 'wrist_left_Z_kin'}));
        else
            ch_wrist = find(ismember({EEG.chanlocs.labels}, ...
                {'wrist_right_X_kin', 'wrist_right_Y_kin', 'wrist_right_Z_kin'}));
        end

        wrist_data = EEG.data(ch_wrist, :);
        Fs = EEG.srate;

        cutoff = 8;                % Hz
        sgolay_window = 31;        % samples
        sgolay_order = 3;          % order 3

        % Signal processing adapted from Vissani 2021

        % === 1. Low-pass filter at 8 Hz ===
        [b,a] = butter(4, cutoff/(Fs/2), 'low');
        wrist_filt = filtfilt(b, a, wrist_data')';      % transpose for filtfilt

        % === 2. Savitzky–Golay smoothing ===
        wrist_smooth = sgolayfilt(double(wrist_filt'), sgolay_order, sgolay_window)'; 

        % === 3. Forward difference ===
        vel_xyz = diff(wrist_smooth, 1, 2) * Fs;  
        vel_abs_raw = sqrt(sum(vel_xyz.^2, 1));
        
        % Moving mean smoothing
        vel_abs = movmean(vel_abs_raw, round(Fs * 0.2)); 
        
        % Elongate signal to maintain original length
        vel_smooth = [vel_abs(1), vel_abs];

        % Initial baseline threshold
        start_sample = EEG.event(1).latency + start_offset;
        baseline_window_s = 500; % ms
        k_std = 4;

        baseline_samples = find(time >= start_sample*2.5 & ...
                                time <= start_sample*2.5 + baseline_window_s);
        mean_baseline = mean(vel_smooth(baseline_samples));
        std_baseline = std(vel_smooth(baseline_samples));
        threshold = mean_baseline + k_std * std_baseline;
        

        % Initialize structure to store events
        events = struct('A', [], 'B', [], 'C', [], 'D', [], 'E', [], 'F', []);
        events_times = struct('A', [], 'B', [], 'C', [], 'D', [], 'E', [], 'F', []);
        
        % Save preprocessed wrist velocity
        kinematic_block = struct();
        kinematic_block.velocity = vel_smooth;

        search_start = start_sample;
        
        % Definition of number of trials per subject
        if any(strcmp(subject_id, {'wue09'})) || ...
           (any(strcmp(subject_id, {'wue10', 'wue05', 'wue03'})) && f == 2) || ...
           (any(strcmp(subject_id, {'wue11', 'wue05'})) && f == 1) || ...
           (any(strcmp(subject_id, {'wue07'})) && (f == 1 || f == 3)) 
            num_trials = 11;
        elseif any(strcmp(subject_id, {'wue10'})) && f == 1
            num_trials = 9;
        elseif any(strcmp(subject_id, {'wue06'})) && f == 1 || ...
               (any(strcmp(subject_id, {'wue11'})) && f == 3)
            num_trials = 12;
        elseif any(strcmp(subject_id, {'wue07', 'wue11'})) && f == 2 || ...
               (any(strcmp(subject_id, {'wue05'})) && f == 3)
            num_trials = 13;
        else
            num_trials = 10;
        end
        
        % Parameters shared across subjects
        A_add = 0.04;
        c_threshold = 0.053;
        
        c_lag_1 = 10;
        c_lag_2 = 30;
        
        E_lag = 600;
        F_add = 0.03;

        for trial = 1:num_trials

            % Event A --> baseline threshold + A_add offset
            A_idx_rel = find(vel_smooth(search_start:end) > threshold + A_add, 1, 'first');
            if isempty(A_idx_rel), break; end
            A_idx = search_start + A_idx_rel - 1;
            events(trial).A = A_idx;

            % Event B
            max_reach_idx = A_idx + 500; % Search bound for B peak
            [~, rel_B] = max(vel_smooth(A_idx:max_reach_idx));
            B_idx = A_idx + rel_B - 1;
            events(trial).B = B_idx;

            % Event C --> 700-sample search bound, c_threshold
            % close to 0.05 (Vissani 2021)
            C_idx_rel = find(vel_smooth(B_idx:B_idx+700) < c_threshold, 1, 'first');
            C_idx = B_idx + C_idx_rel - 1;
            events(trial).C = C_idx;

            % Local baseline (between c_lag_1 and c_lag_2) and D onset
            % (D_add tuned for subjects with higher noise)
            local_base_start = C_idx + c_lag_1;
            local_base_end = C_idx + c_lag_2;
            mean_local = mean(vel_smooth(local_base_start:local_base_end));
            std_local = std(vel_smooth(local_base_start:local_base_end));
            threshold_local = mean_local + k_std * std_local;

            D_idx_rel = find(vel_smooth(C_idx+c_lag_1:end) > threshold_local + D_add, 1, 'first');
            if isempty(D_idx_rel), break; end
            D_idx = C_idx + c_lag_1 + D_idx_rel - 1;
            events(trial).D = D_idx;

            % Event E --> peak within 300 samples from D
            [~, rel_E] = max(vel_smooth(D_idx:D_idx+300));
            E_idx = D_idx + rel_E - 1;
            events(trial).E = E_idx;

            % Event F
            F_idx_rel = find(vel_smooth(E_idx:E_idx+E_lag) < (mean_baseline + F_add), 1, 'first');
            F_idx = E_idx + F_idx_rel - 1;
            events(trial).F = F_idx;

            % Update search starting point
            search_start = F_idx + F_lag;
        end
        
        % === Plot ===
        figure;
        plot(time, vel_smooth, 'k', 'DisplayName', 'Absolute Velocity'); 
        ylim([-1,1]);
        xlabel('Time (ms)');
        ylabel('Absolute velocity (m/s)');
        title(sprintf('Absolute velocity with kinematic events superimposed - Subject %s - Block %d', subject_list{s}, f));
        grid on;
        hold on;
        
        event_names = fieldnames(events);
        colors = lines(length(event_names));
        
        for e = 1:length(event_names)
            ev_name = event_names{e};
            all_idx = []; % collect all indices for this event
        
            for trial = 1:length(events)
                all_idx = [all_idx, events(trial).(ev_name)];
            end
        
            % Single plot per event type
            plot(time(all_idx), vel_smooth(all_idx), 'o', ...
                 'Color', colors(e,:), 'MarkerSize', 8, 'DisplayName', ev_name);
        end
        
        legend('Location','bestoutside');

        % Save high-resolution PNG and .FIG
        set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]);
        set(gca, 'FontSize', 14);
        fig_name = sprintf('%s_block%d', subject_id, f);
        print(fullfile(set_path, fig_name), '-dpng', '-r300'); % -r300 = 300 dpi
        savefig(gcf, fullfile(set_path, [fig_name '.fig']));

        % Build structure with event times in ms (ultimately not used)
        for t = 1:length(events)
            fn = fieldnames(events);
            for k = 1:length(fn)
                if ~isempty(events(t).(fn{k}))
                    events_times(t).(fn{k}) = events(t).(fn{k}) * 2.5; 
                end
            end
        end

    end
end
