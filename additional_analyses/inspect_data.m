% inspect_data.m  -- scratch: understand data structures before building analysis
base = 'F:\Projects\Parkinson_ReachGrasp\Reprocessing';
ep   = fullfile(base,'RESULTS_final','Epochs','Epochs_allSubjects.mat');
S    = load(ep);
fn   = fieldnames(S);
fprintf('TOP VARS: %s\n', strjoin(fn,', '));
E    = S.(fn{1});
fprintf('EPOCHS fields: %s\n', strjoin(fieldnames(E),', '));
disp('--- params ---'); disp(E.params);
subj = E.params.subjects{1};
fprintf('First subj: %s\n', subj);
disp('--- subj fields ---'); disp(fieldnames(E.(subj)));
fprintf('fs=%g\n', E.(subj).fs);
fprintf('onset.data size: %s, n=%d\n', mat2str(size(E.(subj).onset.data)), E.(subj).onset.n);
fprintf('pull.data  size: %s, n=%d\n', mat2str(size(E.(subj).pull.data)),  E.(subj).pull.n);
fprintf('onset.times: [%g %g] npts=%d\n', E.(subj).onset.times(1), E.(subj).onset.times(end), numel(E.(subj).onset.times));
fprintf('pull.times : [%g %g] npts=%d\n', E.(subj).pull.times(1),  E.(subj).pull.times(end),  numel(E.(subj).pull.times));
disp('--- full subj struct ---'); disp(E.(subj));

fprintf('\nwhich eeglab: %s\n', which('eeglab'));

if ~isempty(which('eeglab'))
    eeglab nogui;
    lfpdir = fullfile(base,'wue02','Preprocessed','LFP');
    ff = dir(fullfile(lfpdir,'*_wEv.set'));
    L = pop_loadset('filename', ff(1).name, 'filepath', lfpdir);
    fprintf('\nLFP set: nbchan=%d srate=%g pnts=%d\n', L.nbchan, L.srate, L.pnts);
    fprintf('LFP chan labels: %s\n', strjoin({L.chanlocs.labels}, ', '));

    eegdir = fullfile(base,'wue02','Preprocessed','EEG');
    fe = dir(fullfile(eegdir,'*_manual.set'));
    EEG = pop_loadset('filename', fe(1).name, 'filepath', eegdir);
    fprintf('\nEEG set: nbchan=%d srate=%g\n', EEG.nbchan, EEG.srate);
    fprintf('EEG chan labels: %s\n', strjoin({EEG.chanlocs.labels}, ', '));

    emgdir = fullfile(base,'wue02','01_Extracted','EMG_KIN');
    fm = dir(fullfile(emgdir,'*.set'));
    EMG = pop_loadset('filename', fm(1).name, 'filepath', emgdir);
    fprintf('\nEMG set: nbchan=%d srate=%g\n', EMG.nbchan, EMG.srate);
    fprintf('EMG chan labels: %s\n', strjoin({EMG.chanlocs.labels}, ', '));
end
fprintf('\nINSPECT DONE\n');
