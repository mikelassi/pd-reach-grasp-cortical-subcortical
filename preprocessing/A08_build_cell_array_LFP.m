%% ======================================
% Script: Extract LFP segments, trials, and phases
%
% Purpose:
%   - Extract LFP segments based on task events (A_T[x], A→B, B→C, ..., F→mid(F,A_T[x+1])).
%   - Construct full LFP trials from start to mid(F, A_Tx).
%   - Organize LFP data into phases:
%       • rest_bef
%       • reaching
%       • grasping
%       • pulling
%       • rest_aft
%   - Visualize LFP phases consistent with Vissani 2021: Rest, Reaching, Grasping, Pulling.
%
% Special handling:
%   - Trials of same length enforced for subjects 'wue05', 'wue07', 'wue09', 'wue11' 
%     (from 2.8 s before to 0.7 s after F event – note: currently not used).
%   - Excluded trials are handled explicitly per subject and block.
%   - Subject wue06, block 2 uses a baseline before TENStrigger artifact instead of 50 samples after TENStrigger.
%   - Kinematic data is aligned to LFP trials and adjusted for excluded trials.
%
% Outputs:
%   - LFP_segments: cell array [trial × region]
%   - LFP_phases: cell array [trial × phase]
%   - LFP_trials: concatenated trials [trial × time]
%   - kinematic_4LFP_trials: aligned kinematic velocity per trial
%   - baseline_LFP: baseline segment from start → A_T1 (or TENStrigger region)
%   - trial durations (samples and ms)
%
% Author: Tommaso Marcantoni
% Date: 2026
%% ======================================

clear; clc;

subject_list = {'wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11'};

kin_subfolder = '02_Kinematics\Events';

% Struct for baselines and trial durations
all_baselines = struct();
all_durations = {};

% === Trials rimossi nel preprocessing LFP ===
removed_trials = struct( ...
    'wue11', [1, 5], ...   % file 1, trial 5
    'wue05', [3, 11] ...   % file 3, trial 11
);

[ALLLFP, LFP, CURRENTSET] = eeglab;
ALLLFP = [];

for s = 1:length(subject_list)
    
    subject_id = subject_list{s};
    base_path  = 'H:\Parkinson_ReachGrasp\Reprocessing';
    preproc_path = fullfile(base_path, subject_id, 'Preprocessed', 'LFP');

    all_baselines.(subject_id) = {}; 
    
    original_LFP_path = fullfile(base_path,subject_id, '03_SyncRaw', 'LFP_wEv');
    or_set_files = dir(fullfile(original_LFP_path, '*.set'));
    kin_path = fullfile(base_path, subject_id, kin_subfolder);
    
    set_files = dir(fullfile(preproc_path, '*.set'));
    kin_files = dir(fullfile(kin_path, '*_kinematic_block.mat'));
    
    for f = 1:length(set_files)

        set_file = set_files(f).name;
        fprintf('\nLoad the preprocessed file: %s\n', set_file);
        
        LFP = pop_loadset('filename', set_file, 'filepath', preproc_path);
        LFP_or = pop_loadset('filename', or_set_files(f).name, 'filepath', original_LFP_path);
        Fs_LFP = LFP.srate;  

        event_types = {LFP.event.type};
        event_latencies = round([LFP.event.latency]); % sample indices
        
        event_types_or = {LFP_or.event.type};

        % Determine number of trials based on F_T events 
        is_FT = startsWith(event_types, 'F_T');                                      
        nums = cellfun(@(x) str2double(extractAfter(x, 'F_T')), event_types(is_FT)); 
        [num_trials, ~] = max(nums); % total number of trials

        % === Gestione trial esclusi ===
        if isfield(removed_trials, subject_id) && numel(removed_trials.(subject_id)) >= 2 ...
                && removed_trials.(subject_id)(1) == f
            trial_to_remove = removed_trials.(subject_id)(2);
            fprintf('⚠️  %s file %d: escluso trial %d\n', subject_id, f, trial_to_remove);
        else
            trial_to_remove = [];
        end
        
        % === Lista dei trial validi ===
        valid_trials = setdiff(1:num_trials, trial_to_remove);
        num_valid_trials = numel(valid_trials);
        fprintf('→ Trials validi: %s\n', mat2str(valid_trials));

        % Nuove regioni (7 in totale)
        region_labels = {'start_A','A_B','B_C','C_D','D_E','E_F','F_midFA'};
        num_regions   = length(region_labels);
        
        % Prealloca cell array: trials validi × regions
        LFP_segments = cell(num_valid_trials, num_regions);
        baseline_LFP = []; % Store baseline from beginning to A_T1

        %% === COSTRUZIONE SEGMENTI ===
        for new_t = 1:num_valid_trials
            t = valid_trials(new_t);  % trial originale

            % --- Trova evento A_Tx ---
            A_idx = find(strcmp(event_types, sprintf('A_T%d', t)));
            if isempty(A_idx)
                next_valid = find(~cellfun('isempty', regexp(event_types, 'A_T\d+')));
                next_valid_nums = cellfun(@(x) str2double(extractAfter(x, 'A_T')), event_types(next_valid));
                next_t = min(next_valid_nums(next_valid_nums > t));
                if isempty(next_t)
                    warning('Trial %d non trovato e nessun successivo valido.', t);
                    continue;
                else
                    fprintf('Trial %d mancante → salto al trial %d\n', t, next_t);
                    continue;
                end
            end

            % ---- 1. start → A_T1 ----
            A1_idx = find(strcmp(event_types, sprintf('A_T%d', t)));
            if isempty(A1_idx)
                warning('Trial %d: A_T%d not found.', t, t);
                continue;
            end

            end_latency = event_latencies(A1_idx);

            if t == 1 
                start_latency = 1;
                baseline_LFP = LFP.data(:, start_latency:end_latency);
                all_baselines.(subject_id){f} = baseline_LFP;
            else 
                start_latency = mid_FA + 1;
            end
            
            LFP_segments{new_t,1} = LFP.data(:, start_latency:end_latency);

            % ---- 2-6. regioni canoniche A→B ... E→F ----
            base_regions = {'A','B','C','D','E','F'};
            for r = 1:(length(base_regions)-1)
                start_label = sprintf('%s_T%d', base_regions{r}, t);
                end_label   = sprintf('%s_T%d', base_regions{r+1}, t);
                
                start_idx = find(strcmp(event_types, start_label));
                end_idx   = find(strcmp(event_types, end_label));

                if isempty(start_idx) || isempty(end_idx)
                    warning('Trial %d: %s→%s not found.', t, base_regions{r}, base_regions{r+1});
                    continue;
                end
                LFP_segments{new_t,r+1} = LFP.data(:, event_latencies(start_idx)+1:event_latencies(end_idx));
            end
            
            % ---- 7. F → mid(F, A_T1_next) ----
            F_idx = find(strcmp(event_types, sprintf('F_T%d', t)));

            if ~isempty(F_idx)
                F_lat = event_latencies(F_idx);
            
                % Trova A del trial successivo (se esiste)
                nextA_idx = find(strcmp(event_types, sprintf('A_T%d', t+1)));
            
                if isempty(nextA_idx)
                    % Se il trial successivo è stato rimosso → cerca un "boundary" dopo F
                    boundary_idx = find(strcmp(event_types, 'boundary'));
                    boundary_lat = event_latencies(boundary_idx);
                    boundary_after_F = boundary_lat(boundary_lat > F_lat);
            
                    if ~isempty(boundary_after_F)
                        end_lat = round(boundary_after_F(1));

                        % Verifico se questo boundary è l'ultimo evento del dataset
                        is_last_boundary = end_lat == round(LFP.event(end).latency);

                        end_lat = min(end_lat, size(LFP.data,2));
                    
                        if is_last_boundary
                            % Ultimo trial → chiudo regolarmente al boundary finale
                            mid_FA = end_lat;
                            fprintf('✅ Trial %d: ultimo trial, fine al boundary finale (%.1f s)\n', ...
                                    t, end_lat / Fs_LFP);
                        else
                            % Trial intermedio → A_T mancante, uso boundary come sostituto
                            mid_FA = end_lat;
                            fprintf('⚠️ Trial %d: A_T%d mancante, uso boundary a %.1f s come fine rest_aft\n', ...
                                    t, t+1, end_lat / Fs_LFP);
                        end
                    else
                        % Nessun boundary → fino alla fine del segnale
                        end_lat = LFP.pnts;
                        fprintf('⚠️ Trial %d: A_T%d e boundary mancanti → fino a fine recording\n', ...
                                t, t+1);
                    end
                else
                    % Se A_next esiste → prendo la metà tra F e A_next
                    A_next_lat = event_latencies(nextA_idx);
                    end_lat = round(F_lat + (A_next_lat - F_lat)/2);
                    mid_FA = end_lat;
                end
            
                LFP_segments{new_t,7} = LFP.data(:, F_lat+1:end_lat);
            
            else
                warning('Trial %d: F_T%d non trovato.', t, t);
            end
        end

        fprintf('✅ LFP segments creati: %d trials × %d regions\n', size(LFP_segments,1), size(LFP_segments,2));

        %% === FASI ===
        phase_names   = {'rest_bef','reaching','grasping','pulling', 'rest_aft'};
        phase_regions = {1, [2 3], 4, [5 6], 7};  

        LFP_phases = cell(num_valid_trials, length(phase_names));
        LFP_trials = cell(num_valid_trials, 1);
        trial_durations = zeros(num_valid_trials, 2);

        for new_t = 1:num_valid_trials
            t = valid_trials(new_t);
            A1_idx  = find(strcmp(event_types, sprintf('A_T%d', t)));
            if isempty(A1_idx), continue; end
            A1_latency = event_latencies(A1_idx);
            trial_data = [];

            for p = 1:length(phase_names)
                tmp_data = [];
                for r = phase_regions{p}
                    if new_t == 1 && r == 1
                        start_rest_bef = A1_latency - 250;
                        tmp_data = [tmp_data, LFP_segments{new_t,r}(:, start_rest_bef:end)];
                    elseif r <= size(LFP_segments,2) && ~isempty(LFP_segments{new_t,r})
                        tmp_data = [tmp_data, LFP_segments{new_t,r}];
                    end
                end
                LFP_phases{new_t,p} = tmp_data;
                trial_data = [trial_data, tmp_data];  
            end

            LFP_trials{new_t,1} = trial_data;
            dur_samples = length(trial_data);
            dur_ms = (dur_samples / Fs_LFP) * 1000;

            trial_durations(new_t,:) = [dur_samples, dur_ms];
            all_durations(end+1,:) = {subject_id, f, new_t, dur_samples, dur_ms};
        end

        fprintf('✅ LFP phases create: %d valid trials × %d phases\n', size(LFP_phases,1), size(LFP_phases,2));
        
        %% === Carica kinematica ===
        kin_file = kin_files(f).name;
        
        if exist(fullfile(kin_path, kin_file), 'file')
            load(fullfile(kin_path, kin_file), 'kinematic_block'); % Load smooth velocity for current block
        else
            warning('Nessun file cinematico trovato per %s - %s', subject_id, set_file);
            continue;
        end

    % ========================================================
        % === NEW ROBUST KINEMATIC ALIGNMENT (NO TENS NEEDED) ===
        % ========================================================
        
        % 1. Find A_T in the ORIGINAL file (FORCE FIRST MATCH ONLY)
        A1_idx_or = find(strcmp({LFP_or.event.type}, 'A_T1'), 1, 'first');
        if isempty(A1_idx_or)
            A1_idx_or = find(startsWith({LFP_or.event.type}, 'A_T'), 1, 'first');
        end
        A1_latency_or = round(LFP_or.event(A1_idx_or).latency);
        
        % 2. Find A_T in the PREPROCESSED file (FORCE FIRST MATCH ONLY)
        A1_idx_prep = find(strcmp({LFP.event.type}, 'A_T1'), 1, 'first');
        if isempty(A1_idx_prep)
            A1_idx_prep = find(startsWith({LFP.event.type}, 'A_T'), 1, 'first');
        end
        A1_latency_prep = round(LFP.event(A1_idx_prep).latency);
        
        % 3. Calculate EXACT Kinematic start and end indices using relative offset
        samples_chopped = A1_latency_or - A1_latency_prep;
        
        start_idx_kin = samples_chopped + 1;
        end_idx_kin   = start_idx_kin + size(LFP.data, 2) - 1;
        
        % 4. Safely extract the matching velocity block
        if end_idx_kin <= length(kinematic_block.velocity)
            vel_block = kinematic_block.velocity(start_idx_kin:end_idx_kin);
        else
            warning('Kinematic block is shorter than LFP. Truncating to fit.');
            vel_block = kinematic_block.velocity(start_idx_kin:end);
        end
        
        time = (0:length(vel_block)-1) / Fs_LFP;
        % ========================================================
        % ========================================================
        
        if (strcmp(subject_id, 'wue11') && f == 1)
            trial_to_remove = 5;
        elseif (strcmp(subject_id, 'wue05') && f == 3)
            trial_to_remove = 11;
        else
            trial_to_remove = [];
        end
        
        event_latencies_or = round([LFP_or.event.latency]);
        
        if ~isempty(trial_to_remove)
            % Trova eventi delimitanti la finestra da saltare
            F_idx_pre  = find(strcmp(event_types_or, sprintf('F_T%d', trial_to_remove - 1)));
            A_idx_curr = find(strcmp(event_types_or, sprintf('A_T%d', trial_to_remove)));
            F_idx_curr = find(strcmp(event_types_or, sprintf('F_T%d', trial_to_remove)));
            A_idx_next = find(strcmp(event_types_or, sprintf('A_T%d', trial_to_remove + 1)));
        
            F_lat_pre  = event_latencies_or(F_idx_pre);
            A_lat_curr = event_latencies_or(A_idx_curr);
            F_lat_curr = event_latencies_or(F_idx_curr);
            A_lat_next = event_latencies_or(A_idx_next);
        
            mid_FA_pre  = round(F_lat_pre + (A_lat_curr - F_lat_pre)/2);
            mid_FA_post = round(F_lat_curr + (A_lat_next - F_lat_curr)/2);
        
            % Converti in indici cinematici (usando il nuovo offset calcolato)
            kin_start = mid_FA_pre - samples_chopped;   % prima di mid_FA
            kin_end   = mid_FA_post - samples_chopped;  % fino a mid_FA post
        
        
            vel_block(kin_start:kin_end) = [];
        end
        
        % === Estrarre trial validi ===
        kinematic_4LFP_trials = cell(num_valid_trials, 1);
        offset = A1_latency_prep - 250;
        start_idx = offset;
        
        for new_t = 1:num_valid_trials
            t = valid_trials(new_t);
            trial_len = size(LFP_trials{new_t}, 2);
        
            end_idx = start_idx + trial_len - 1;
        
            if end_idx <= length(vel_block)
                trial_vel = vel_block(start_idx:end_idx);
                kinematic_4LFP_trials{new_t} = trial_vel;
            else
                warning('⚠️ Trial %d (%s): cinematico troppo corto (end_idx=%d, vel_len=%d)', ...
                    t, subject_id, end_idx, length(vel_block));
                kinematic_4LFP_trials{new_t} = vel_block(start_idx:end);
            end
        
            start_idx = end_idx + 1;
        end


        %% === SAVE ===
        save(fullfile(preproc_path, [set_file(1:end-4) '_LFP_trialsByRegionAndPhase.mat']), ...
            'LFP_trials','LFP_segments','LFP_phases','region_labels','phase_names','phase_regions','kinematic_4LFP_trials', 'baseline_LFP');
    end
end

% % === File TXT globale per durata trials ===
% out_file = fullfile(base_path, 'RESULTS/OLD/Duration_trials&baseline_01/', 'all_trials_durations.txt');
% fid = fopen(out_file,'w');
% fprintf(fid, 'Subject\tBlock\tTrial\tDur_samples\tDur_ms\n');
% for i = 1:size(all_durations,1)
%     fprintf(fid, '%s\t%d\t%d\t%d\t%.2f\n', ...
%         all_durations{i,1}, all_durations{i,2}, all_durations{i,3}, ...
%         all_durations{i,4}, all_durations{i,5});
% end
% fclose(fid);
