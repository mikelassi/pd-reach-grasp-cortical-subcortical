%% LFP PREPROCESSING PIPELINE (EEGLAB)
% ------------------------------------------------------------
% Project: Cortical–Subcortical Connectivity in Parkinson’s Disease
% Author: Tommaso Marcantoni
%
% Description:
%   Preprocessing pipeline for STN LFP recordings during a
%   reach-and-grasp task.
%
%   The pipeline includes:
%     - Event-based temporal cropping
%     - Subject-specific trial rejection
%     - Cardiac (QRS) artifact template subtraction (subject wue03)
%     - Channel-wise z-score normalization
%
% Dependencies:
%   - EEGLAB (tested with v2025.0.0)
%   - Signal Processing Toolbox
%
% Input:
%   - EEGLAB .set files with events:
%       TENStrigger, A_Tx, F_Tx
%
% Output:
%   - Preprocessed LFP datasets (.set)
%
% Notes:
%   - Dataset-specific preprocessing
%   - Manual and subject-dependent corrections are applied
%   - Script provided as raw reference implementation
% ------------------------------------------------------------

clear; clc;

subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'};
%subject_list = {'wue03'};


subject_info = {};
[ALLLFP, LFP, CURRENTSET] = eeglab;
ALLLFP = [];

%% PATHS (EDIT THIS SECTION)

base_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';

% IMPORTANT:
% Subject wue03 presents a strong cardiac (QRS) artifact in LFP channel 2.
% A subject-specific preprocessing branch is applied:
%   - QRS peak detection
%   - Epoching around QRS
%   - Global artifact template estimation
%   - Template subtraction from the LFP signal

for s = 1:length(subject_list)

    subject_id = subject_list{s};
    set_path = fullfile(base_path, subject_id, 'Extracted', 'LFP_wEV'); % Path dei .set

    % Lista dei file .set
    set_files = dir(fullfile(set_path, '*.set'));

    % If we consider subject wue03 we need to find QRS complex,
    % epoch data and subtract mean artefact to the signal
    if any(strcmp(subject_id, {'wue03'})) 

        eeg_set_path = fullfile(base_path, subject_id, 'Extracted', 'EEG_wEV'); 
        eeg_set_files = dir(fullfile(eeg_set_path, '*.set'));

        all_epochs_allblocks = []; 
        all_locs = cell(1, length(set_files));

        for f = 1:length(set_files)

            set_file = set_files(f).name;
            fprintf('\nCarico file LFP: %s\n', set_file);
            
            LFP = pop_loadset('filename', set_file, 'filepath', set_path);
            [~, name_base] = fileparts(set_file);

            % Manual division of signal between 50 sample after TENSTRIGGER  and FT10 + 300 samples
    
            % All events 
            event_types = {LFP.event.type};
    
            % Latenza di TENSTRIGGER
            lat_TEN_idx_all = find(strcmp(event_types, 'TENStrigger'));
    
            % Latenza A_T1
            AT1_latency_idx = find(strcmp(event_types, 'A_T1'));
            AT1_latency = LFP.event(AT1_latency_idx).latency;
            
            % First TENSTRIGGER
            lat_TEN_idx = lat_TEN_idx_all(1);
            lat_TEN = LFP.event(lat_TEN_idx).latency;
    
            % Find event F_Tx latency (first need to find x--> n° of trials)
            is_FT = startsWith(event_types, 'F_T');                                      % all events that starts with F_T
            nums = cellfun(@(x) str2double(extractAfter(x, 'F_T')), event_types(is_FT)); % extract num after F_T
            [max_num, idx_max] = max(nums);                                              % find maximum trial (eg. 10)
            last_event_label = ['F_T' num2str(max_num)];
            second_latency_idx = find(strcmp(event_types, ['F_T' num2str(max_num)]));
            second_latency = LFP.event(second_latency_idx).latency;                      % finally F_Tx latency
            
            start_point = lat_TEN + 50;  % TENStrigger latency + 50 samples
            end_point = second_latency + 300;

            % Store baseline informations
            baseline_samples = AT1_latency - start_point;              
            baseline_ms      = (baseline_samples / LFP.srate) * 1000;  
            subject_info(end+1,:) = {subject_id, set_file, baseline_samples, baseline_ms};

            % === Select the region (pop_select) ===  This before resampling to
            % assure same time duration wrt EEG
            original_LFP = LFP;
            LFP = pop_select(LFP, 'point', [(start_point) (end_point)]);

            %% ==============================
            % CARDIAC (QRS) ARTIFACT ESTIMATION
            % ==============================

            % QRS detection --> simple findpeaks after band pass filtering 
            lowcut  = 3;    % Hz  → elimina drift lenti
            highcut = 50;   % Hz  → elimina rumore ad alta frequenza
            [b,a] = butter(4, [lowcut highcut] / (LFP.srate/2), 'bandpass');

            LFP_filt = filtfilt(b,a, double(LFP.data(2,:)));

            min_threshold = mean(LFP_filt)+2*std(LFP_filt);
            
            % Find location of negative peaks
            search_start = 100; % offset TeNStrigger signal
            signal_to_search = LFP_filt(search_start:end);
            [peaks,locs] = findpeaks(-signal_to_search, "MinPeakDistance",round(0.5*LFP.srate),"MinPeakHeight",-min_threshold);
            peaks = -peaks; % Conversione a valori originali
            locs = locs + search_start;
            
            % Plot of LFP signal with peaks superimposed
            % times
            t = LFP.times;

            % Plot LFP
            figure;
            plot(t, LFP.data(2,:), 'b'); hold on;

            % Plot R-peaks
            plot(locs*2.5, LFP.data(2,locs), 'ro', 'MarkerFaceColor','r');

            xlabel('Time (ms)');
            ylabel('Amplitude (mV)');
            title('ECG with detected QRS peaks');
            legend('ECG','R-peaks');
        
            % Epoching around the peak
            epoch_ms = [-100 200];   % -200 ms to +300 ms
            samples_per_ms = LFP.srate / 1000;   % samples per ms

            epoch_pts = round(epoch_ms * samples_per_ms);   % Interval in samples

            all_epochs = [];  

            % for loop on number of peaks
            for k = 1:length(locs)

                % indexes of the start and end of epochs
                idx_start = locs(k) + epoch_pts(1);
                idx_end   = locs(k) + epoch_pts(2);

                if idx_start < 1 || idx_end > size(LFP.data,2)
                    continue;  
                end

                % Extract signal in the epoch
                epoch = LFP.data(2, idx_start:idx_end);
                all_epochs = [all_epochs; epoch]; 
            end

            all_epochs_allblocks = [all_epochs_allblocks; all_epochs];  
            all_locs{f} = locs;

        end

        for f = 1:length(set_files)

            set_file = set_files(f).name;
            fprintf('\nCarico file LFP: %s\n', set_file);

            LFP = pop_loadset('filename', set_file, 'filepath', set_path);
            [~, name_base] = fileparts(set_file);

            % Manual division of signal between 50 sample after TENSTRIGGER  and FT10 + 300 samples

            % All events 
            event_types = {LFP.event.type};

            % Latency of TENSTRIGGER
            lat_TEN_idx_all = find(strcmp(event_types, 'TENStrigger'));

            % Latency of A_T1
            AT1_latency_idx = find(strcmp(event_types, 'A_T1'));
            AT1_latency = LFP.event(AT1_latency_idx).latency;

            % First TENSTRIGGER
            lat_TEN_idx = lat_TEN_idx_all(1);
            lat_TEN = LFP.event(lat_TEN_idx).latency;

            % Find event F_Tx latency (first need to find x--> n° of trials)
            is_FT = startsWith(event_types, 'F_T');                                      % all events that starts with F_T
            nums = cellfun(@(x) str2double(extractAfter(x, 'F_T')), event_types(is_FT)); % extract num after F_T
            [max_num, idx_max] = max(nums);                                              % find maximum trial (eg. 10)
            last_event_label = ['F_T' num2str(max_num)];
            second_latency_idx = find(strcmp(event_types, ['F_T' num2str(max_num)]));
            second_latency = LFP.event(second_latency_idx).latency;                      % finally F_Tx latency

            start_point = lat_TEN + 50;  % TENStrigger latency + 50 samples
            end_point = second_latency + 300;

            % === Select the region (pop_select) ===  This before resampling to
            % assure same time duration wrt EEG
            original_LFP = LFP;
            LFP = pop_select(LFP, 'point', [(start_point) (end_point)]);

            global_mean_epoch = mean(all_epochs_allblocks,1);

            epoch_time = linspace(epoch_ms(1), epoch_ms(2), size(all_epochs,2));
            
            % Just to plot one time only
            if f == 1
                figure; hold on;

                for k = 1:size(all_epochs_allblocks,1)
                    plot(epoch_time, all_epochs_allblocks(k,:), 'b');
                end


                plot(epoch_time, global_mean_epoch, 'k', 'LineWidth',2);

                xlabel('Time (ms)');
                ylabel('Amplitude');
                title('Global mean artifact across all blocks');

                % Save in PNG HD
                save_path = fullfile(base_path, subject_id, 'Preprocessed_02', 'LFP_01');
                if ~exist(save_path, 'dir')
                    mkdir(save_path);
                end
                set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1200, 600]); 
                set(gca, 'FontSize', 12); 
                fig_name = sprintf('template_LFP_sub_%s_block%d', subject_id, f);
                print(fullfile(save_path, fig_name), '-dpng', '-r300');
            end


            % Creation of the artefactual signal to be subtracted
            artifact_signal = zeros(1, size(LFP.data,2));
            locs = all_locs{f};

            for k = 1:length(locs)
                idx_start = locs(k) + epoch_pts(1);
                idx_end   = locs(k) + epoch_pts(2);

                if idx_start < 1 || idx_end > size(LFP.data,2)
                    continue;
                end

                artifact_signal(idx_start:idx_end) = artifact_signal(idx_start:idx_end) + global_mean_epoch;
            end

            %% ==============================
            % CARDIAC ARTIFACT TEMPLATE SUBTRACTION
            % ==============================

            % Subtract the artifactual signal
            LFP_corrected = LFP.data(2,:) - artifact_signal;

            % Visualization of original and corrected signal
            t_start = 30; % sec
            t_end   = 35; % sec

            % Times in s
            t_sec = LFP.times / 1000;
            idx_win = t_sec >= t_start & t_sec <= t_end;

            figure;
            plot(t_sec(idx_win), LFP.data(2, idx_win), 'b'); hold on;
            plot(t_sec(idx_win), LFP_corrected(idx_win), 'r');
            legend('Original LFP','Corrected LFP');
            xlabel('Time (s)');
            ylabel('Amplitude');
            title(sprintf('Cardiac artifact removal (%.1f–%.1f s)', t_start, t_end));

            % Update the final LFP structure with both LFP channels (one
            % only modified, the second one)
            LFP.data(2,:) = LFP_corrected;

            if (strcmp(subject_id, 'wue03') && f == 3)
                trial_to_remove = 4;
            else
                trial_to_remove = [];
            end

            if ~isempty(trial_to_remove)
            
                % --- Extract latencies and events ---
                event_types = {LFP.event.type};
                event_latencies = round([LFP.event.latency]);
            
                % Beginning and end latencies of trials to remove 
                F_idx_pre  = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove - 1)));
                A_idx_curr = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove)));
                F_idx_curr = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove)));
                A_idx_next = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove + 1)));
            
                % Calculation of mean points to define window to remove
                F_lat_pre  = event_latencies(F_idx_pre);
                A_lat_curr = event_latencies(A_idx_curr);
                F_lat_curr = event_latencies(F_idx_curr);
                A_lat_next = event_latencies(A_idx_next);
        
                mid_FA_pre  = round(F_lat_pre + (A_lat_curr - F_lat_pre)/2) - 1;
                mid_FA_post = round(F_lat_curr + (A_lat_next - F_lat_curr)/2);
        
                % --- Window removal ---
                LFP = pop_select(LFP, 'rmpoint', [mid_FA_pre mid_FA_post]);
            end

            %% ==============================
            % NORMALIZATION AND SAVE
            % ==============================

            % === Z-score normalization per channel ===
            LFP.data = zscore(LFP.data, 0, 2);

            % === Save preprocessed dataset ===
            out_path = fullfile(base_path, subject_id, 'Preprocessed_02', 'LFP_01');
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            LFP = pop_saveset(LFP, 'filename', ['preproc_' set_file], 'filepath', out_path);

            % Store in memory EEGLAB
            [ALLLFP, LFP, CURRENTSET] = eeg_store(ALLLFP, LFP, f);
        end
    else

    % NOTE:
    % For all other subjects, no explicit cardiac artifact correction
    % was required based on visual inspection.
    
        for f = 1:length(set_files)
            
            set_file = set_files(f).name;
            fprintf('\nCarico file LFP: %s\n', set_file);
            
            LFP = pop_loadset('filename', set_file, 'filepath', set_path);
            [~, name_base] = fileparts(set_file);
    
            % Manual division of signal between 50 sample after TENSTRIGGER  and FT10 + 300 samples
    
            % All events 
            event_types = {LFP.event.type};
    
            % Latency of TENSTRIGGER
            lat_TEN_idx_all = find(strcmp(event_types, 'TENStrigger'));
    
            % Latency A_T1
            AT1_latency_idx = find(strcmp(event_types, 'A_T1'));
            AT1_latency = LFP.event(AT1_latency_idx).latency;
            
            % First TENSTRIGGER
            lat_TEN_idx = lat_TEN_idx_all(1);
            lat_TEN = LFP.event(lat_TEN_idx).latency;
    
            % Find event F_Tx latency (first need to find x--> n° of trials)
            is_FT = startsWith(event_types, 'F_T');                                      % all events that starts with F_T
            nums = cellfun(@(x) str2double(extractAfter(x, 'F_T')), event_types(is_FT)); % extract num after F_T
            [max_num, idx_max] = max(nums);                                              % find maximum trial (eg. 10)
            last_event_label = ['F_T' num2str(max_num)];
            second_latency_idx = find(strcmp(event_types, ['F_T' num2str(max_num)]));
            second_latency = LFP.event(second_latency_idx).latency;                      % finally F_Tx latency
            
            start_point = lat_TEN + 50;  % TENStrigger latency + 50 samples
            end_point = second_latency + 300;
    
            % Store baseline info
            baseline_samples = AT1_latency - start_point;              
            baseline_ms      = (baseline_samples / LFP.srate) * 1000;  
            subject_info(end+1,:) = {subject_id, set_file, baseline_samples, baseline_ms};
            
    
            % === Select the region (pop_select) ===  This before resampling to
            % assure same time duration wrt EEG
            original_LFP = LFP;
            LFP = pop_select(LFP, 'point', [(start_point) (end_point)]);
            
            event_types = {LFP.event.type};
            event_latencies = round([LFP.event.latency]); % sample indices
    
            % Handle corrupted trials for wue11 ( block 1 trial 5) and wue05 (block 3 trial 11) 
            if (strcmp(subject_id, 'wue11') && f == 1)
                trial_to_remove = 5;
            elseif (strcmp(subject_id, 'wue05') && f == 3)
                trial_to_remove = 11;
            elseif (strcmp(subject_id, 'wue05') && f == 1)
                trial_to_remove = 4;
            else
                trial_to_remove = [];
            end
        
            
            if ~isempty(trial_to_remove)
            
                % --- Extract events and latencies ---
                event_types = {LFP.event.type};
                event_latencies = round([LFP.event.latency]);
            
                % Identify events od start and end to remove
                F_idx_pre  = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove - 1)));
                A_idx_curr = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove)));
                F_idx_curr = find(strcmp(event_types, sprintf('F_T%d', trial_to_remove)));
                A_idx_next = find(strcmp(event_types, sprintf('A_T%d', trial_to_remove + 1)));
            
                % Mean points for removal
                F_lat_pre  = event_latencies(F_idx_pre);
                A_lat_curr = event_latencies(A_idx_curr);
                F_lat_curr = event_latencies(F_idx_curr);
                A_lat_next = event_latencies(A_idx_next);

                if (strcmp(subject_id, 'wue05') && f == 1)
                    mid_FA_pre  = round(F_lat_pre + (A_lat_curr - F_lat_pre)/2) - 1;
                    mid_FA_post = round(F_lat_curr + (A_lat_next - F_lat_curr)/2);
                else
                    mid_FA_pre  = round(F_lat_pre + (A_lat_curr - F_lat_pre)/2);
                    mid_FA_post = round(F_lat_curr + (A_lat_next - F_lat_curr)/2) - 1;
                end
        
                % --- Window removal---
                LFP = pop_select(LFP, 'rmpoint', [mid_FA_pre mid_FA_post]);
            end

            %% ==============================
            % NORMALIZATION AND SAVE
            % ==============================
    
            % === Z-score normalization per channel ===
            LFP.data = zscore(LFP.data, 0, 2);
    
            % === Save dataset preproc ===
            out_path = fullfile(base_path, subject_id, 'Preprocessed_02', 'LFP_01');
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end
            LFP = pop_saveset(LFP, 'filename', ['preproc_' set_file], 'filepath', out_path);
    
            % Store in memory EEGLAB
            [ALLLFP, LFP, CURRENTSET] = eeg_store(ALLLFP, LFP, f);
    
        end
    end
end
