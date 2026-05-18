%% Add events from .mat files to each subject's LFP sets
% The .mat files with event labels (from EMG_KIN) must be copied
% into the same folder as the LFP .set files (in the same order).

clear; clc;

% List of subjects
subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'};

base_path = 'H:\Parkinson_ReachGrasp\Reprocessing\';

for s = 1:length(subject_list)
    
    subject_id = subject_list{s};
    fprintf('\n=== Processing subject: %s ===\n', subject_id);

    set_path = fullfile(base_path, subject_id, '01_Extracted', 'LFP');
    mat_path = fullfile(base_path, subject_id, '02_Kinematics','Events');
    % List .set files
    set_files = dir(fullfile(set_path, '*.set'));

    % List .mat files (with events)
    mat_files = dir(fullfile(mat_path, '*.mat'));

    for f = 1:length(set_files)
        set_file = set_files(f).name;
        mat_file = mat_files(f).name;
        
        fprintf('\nLoading LFP file: %s\n', set_file);
        fprintf('Loading events file: %s\n', mat_file);

        % Load LFP dataset (.set)
        LFP = pop_loadset('filename', set_file, 'filepath', set_path);
        [~, name_base] = fileparts(set_file);

        % Load events (.mat structure)
        events_matrix = load(fullfile(mat_path, mat_file));

        % Number of trials (usually 10)
        n_trials = length(events_matrix.events); 
        event_labels = {'A','B','C','D','E','F'}; 
        n_events = length(event_labels);

        % Convert struct -> numeric matrix [trials x events]
        numeric_events = nan(n_trials, n_events);
        for trial = 1:n_trials
            for e = 1:n_events
                numeric_events(trial, e) = events_matrix.events(trial).(event_labels{e});
            end
        end

        % Build new events structure
        new_events = [];
        for trial = 1:n_trials
            for e = 1:n_events
                ev_latency = numeric_events(trial, e);
                if ~isnan(ev_latency) && ~isempty(ev_latency)
                    tmp_event.type     = sprintf('%s_T%d', event_labels{e}, trial);
                    tmp_event.latency  = ev_latency; % in samples
                    if isempty(LFP.event)
                        tmp_event.urevent = 1;
                    else
                        tmp_event.urevent = LFP.event(end).urevent + length(new_events) + 1;
                    end
                    % Default EEGLAB event fields
                    tmp_event.code     = 'added';
                    tmp_event.bvmknum  = []; 
                    tmp_event.bvtime   = [];
                    tmp_event.channel  = 1;
                    tmp_event.duration = [];
                    tmp_event.visible  = [];
                    new_events = [new_events, tmp_event];
                end
            end
        end

        % Add events to dataset
        LFP.event = [LFP.event, new_events];
        LFP = eeg_checkset(LFP);

        % Save new dataset
        out_path = fullfile(base_path, subject_id, '03_SyncRaw', 'LFP_wEv');
        if ~exist(out_path, 'dir'); mkdir(out_path); end
        out_file = [name_base, '_wEv.set'];
        pop_saveset(LFP, 'filename', out_file, 'filepath', out_path);
        fprintf('Saved: %s\n', out_file);
    end
end
