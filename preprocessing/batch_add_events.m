%% ==============================================================
%  BATCH SCRIPT: Inject Kinematic Events into EEG and LFP datasets
% ===============================================================

clear; clc;

subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'};
base_path  = 'H:\Parkinson_ReachGrasp\Processed\Task\';

% We want to inject events into both modalities to keep them perfectly synced
modalities = {'EEG', 'LFP'};

[ALLEEG, EEG, CURRENTSET] = eeglab;

for s = 1:length(subject_list)
    subject_id = subject_list{s};
    fprintf('\n========================================\n');
    fprintf('=== Injecting Events: Subject %s ===\n', subject_id);
    fprintf('========================================\n');
    
    kin_path = fullfile(base_path, subject_id, 'Preprocessed', 'EMG_KIN');
    
    for m = 1:length(modalities)
        modality = modalities{m};
        
        set_path = fullfile(base_path, subject_id, 'Preprocessed', modality);
        out_path = fullfile(base_path, subject_id, 'Preprocessed', [modality '_wEv']);
        
        if ~exist(out_path, 'dir')
            mkdir(out_path);
        end
        
        set_files = dir(fullfile(set_path, '*.set'));
        
        if isempty(set_files)
            warning('No %s files found for %s. Skipping modality.', modality, subject_id);
            continue;
        end
        
        for f = 1:length(set_files)
            set_file = set_files(f).name;
            [~, name_base] = fileparts(set_file);
            
            % 1. Reconstruct the matching kinematic .mat filename
            % Replaces 'eeg_' or 'lfp_' prefix with 'kinematic_'
            prefix_to_remove = [lower(modality), '_']; 
            kin_base_name = strrep(set_file, prefix_to_remove, 'kinematic_');
            kin_base_name = strrep(kin_base_name, '.set', '');
            
            % Aggressively strip out all known pipeline modifiers
            kin_base_name = strrep(kin_base_name, 'preproc_', '');
            kin_base_name = strrep(kin_base_name, '_wEv', '');
            kin_base_name = strrep(kin_base_name, '_manual', '');
            kin_base_name = strrep(kin_base_name, 'eeg_', '');
            kin_base_name = strrep(kin_base_name, 'lfp_', '');
            
            kin_base_name = strrep(kin_base_name, '_manual', '');

            mat_file = [kin_base_name, '_kinematic_block.mat'];
            
            mat_filepath = fullfile(kin_path, mat_file);
            
            if ~exist(mat_filepath, 'file')
                warning('Missing Kinematic data! Cannot find %s for %s.', mat_file, set_file);
                continue;
            end
            
            fprintf('\n  > Modality: %s | File: %s\n', modality, set_file);
            
            % 2. Load EEG/LFP and Kinematic Events
            EEG = pop_loadset('filename', set_file, 'filepath', set_path);
            events_matrix = load(mat_filepath);
            
            if ~isfield(events_matrix, 'events')
                warning('The file %s does not contain an "events" structure. Skipping.', mat_file);
                continue;
            end
            
            n_trials = length(events_matrix.events);
            event_labels = {'A','B','C','D','E','F'};
            n_events = length(event_labels);
            
            % 3. Extract and Build New Events
            new_events = [];
            event_fields = fieldnames(EEG.event); % Steal the exact structure layout
            
            for trial = 1:n_trials
                for e = 1:n_events
                    ev_latency = events_matrix.events(trial).(event_labels{e});
                    
                    if ~isnan(ev_latency) && ~isempty(ev_latency)
                        
                        % Create a blank event that perfectly matches the BrainVision structure
                        tmp_event = struct();
                        for fld = 1:length(event_fields)
                            tmp_event.(event_fields{fld}) = [];
                        end
                        
                        % Fill in our specific kinematic data
                        tmp_event.type = sprintf('%s_T%d', event_labels{e}, trial);
                        tmp_event.latency = ev_latency; % in samples
                        
                        % Safely add 'code' only if the original file uses it
                        if isfield(tmp_event, 'code')
                            tmp_event.code = 'added_kin'; 
                        end
                        
                        % Add to our new array
                        new_events = [new_events, tmp_event]; 
                    end
                end
            end
            
            % 4. Inject, Check, and Save
            if ~isempty(new_events)
                % Now they match perfectly, MATLAB will happily combine them!
                EEG.event = [EEG.event, new_events]; 
                
                % eeg_checkset with 'makeur' safely rebuilds the urevent table
                EEG = eeg_checkset(EEG, 'makeur'); 
                
                out_file = [name_base, '_wEv.set'];
                pop_saveset(EEG, 'filename', out_file, 'filepath', out_path);
                fprintf('    [OK] Injected %d events. Saved to %s_wEv\n', length(new_events), modality);
            else
                warning('    [FAIL] No valid events found to inject for %s.', set_file);
            end
            % 4. Inject, Check, and Save
            if ~isempty(new_events)
                EEG.event = [EEG.event, new_events];
                % eeg_checkset with 'makeur' safely rebuilds the urevent table
                EEG = eeg_checkset(EEG, 'makeur'); 
                
                out_file = [name_base, '_wEv.set'];
                pop_saveset(EEG, 'filename', out_file, 'filepath', out_path);
                fprintf('    [OK] Injected %d events. Saved to %s_wEv\n', length(new_events), modality);
            else
                warning('    [FAIL] No valid events found to inject for %s.', set_file);
            end
            
            % Clear memory
            ALLEEG = [];
            CURRENTSET = 0;
        end
    end
end
disp('Event injection completed!');