%% Add events from .mat file to subject folder's sets
% Before adding events to set the matlab files containing events labels and
% index have to be dragged from the EMG_KIN folder to EEG (mat files are
% ordered according to set files)

clear; clc;

subject_id = 'wue02';
base_path = 'C:\Users\tomma\OneDrive - University of Pisa\Desktop\TESI\Dataset_tesi';
set_path = fullfile(base_path, subject_id, 'Extracted', 'EEG');

% Set list in folder
set_files = dir(fullfile(set_path, '*.set'));

% Mat list in folder
mat_files = dir(fullfile(set_path, '*.mat'));


[ALLEEG, EEG, CURRENTSET] = eeglab;

for f = 1:length(set_files)
    set_file = set_files(f).name;
    fprintf('\nCarico file: %s\n', set_file);

    mat_file = mat_files(f).name;
    fprintf('\nCarico file: %s\n', mat_file);

    % Load .set file
    EEG = pop_loadset('filename', set_file, 'filepath', set_path);

    [~, name_base] = fileparts(set_file);
    
    %Load the .mat file
    events_matrix = load(mat_file);
    
    n_trials = length(events_matrix.events); % 10 trials normally
    event_labels = {'A','B','C','D','E','F'};
    n_events = length(event_labels);

    % Initialize a matrix to reshape .mat struct for easy handling
    numeric_events = nan(n_trials, n_events);

    for trial = 1:n_trials
        for e = 1:n_events
            numeric_events(trial, e) = events_matrix.events(trial).(event_labels{e});
        end
    end

     
    % Built the struct new_events with 'EEG.event shape'
    new_events = [];
    for trial = 1:n_trials
        for e = 1:n_events
            ev_latency = numeric_events(trial, e);
            if ~isnan(ev_latency) && ~isempty(ev_latency)
                tmp_event.type = sprintf('%s_T%d', event_labels{e}, trial);
                tmp_event.latency = ev_latency; % in samples
                tmp_event.urevent = EEG.event(end).urevent + length(new_events) + 1;
                tmp_event.code = sprintf('added');     %From here below --> default parameters
                tmp_event.bvmknum = []; 
                tmp_event.bvtime =  [];
                tmp_event.channel =  1;
                tmp_event.duration =  [];
                tmp_event.visible =  [];
                new_events = [new_events, tmp_event]; 
            end
        end
    end

    % Add events to EEG struct
    EEG.event = [EEG.event, new_events];
    EEG = eeg_checkset(EEG);

    % Save new set 
    out_file = [name_base, '_wEv.set'];
    out_path = fullfile(base_path, subject_id, 'Extracted', 'EEG_wEv');
    pop_saveset(EEG, 'filename', out_file, 'filepath', out_path);
    fprintf('Salvato: %s\n', out_file);
end

