% diag_emglabels.m -- list EMG channel labels per subject + which my rule selects
base='F:\Projects\Parkinson_ReachGrasp\Reprocessing';
subs={'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
dom={'wue02','wue03','wue10','wue11'};
musc={'iod','triceps','deltoid'};
rp=fullfile(fileparts(which('eeglab')),'plugins','Biosig3.8.4','biosig','maybe-missing'); if exist(rp,'dir'),rmpath(rp);end
eeglab nogui;
for s=1:numel(subs)
    subj=subs{s}; ep=fullfile(base,subj,'01_Extracted','EMG_KIN'); ff=dir(fullfile(ep,'*.set'));
    if isempty(ff), fprintf('%s: no EMG\n',subj); continue; end
    E=pop_loadset('filename',ff(1).name,'filepath',ep);
    lab=lower({E.chanlocs.labels});
    emglab=lab(contains(lab,'_emg'));
    fprintf('\n=== %s  (dominant-rule=%d) ===\n', subj, ismember(subj,dom));
    fprintf('   EMG channels: %s\n', strjoin(emglab,', '));
    if ismember(subj,dom), sel=contains(lab,'_dominant_emg') & ~contains(lab,'non_dominant');
    else, sel=contains(lab,'_emg'); end
    ch_all=find(sel);
    fprintf('   RULE selects: %s\n', strjoin(lab(ch_all),', '));
    for m=1:numel(musc)
        hit=ch_all(contains(lab(ch_all),musc{m}));
        if isempty(hit), fprintf('     %-8s -> NONE\n',musc{m});
        else, fprintf('     %-8s -> %s\n',musc{m}, lab{hit(1)}); end
    end
end
