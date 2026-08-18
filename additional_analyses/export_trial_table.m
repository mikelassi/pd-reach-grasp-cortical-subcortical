clear; clc;
SUBJECTS={'wue02','wue03','wue05','wue06','wue07','wue09','wue10','wue11'};
CONTRA=containers.Map({'wue02','wue03'},{2,2}); BASE='F:\Projects\Parkinson_ReachGrasp\Reprocessing'; fs=400;
rows={}; 
hdr={'subj','block','trial','reach_dur_s','grasp_dur_s','pull_dur_s','mvt_dur_s',...
     'peak_reach_vel','peak_pull_vel','reach_beta_db','grip_beta_db','reach_lowbeta_db','reach_highbeta_db'};
for s=1:numel(SUBJECTS)
  subj=SUBJECTS{s}; if isKey(CONTRA,subj), ch=CONTRA(subj); else, ch=1; end
  pl=fullfile(BASE,subj,'Preprocessed','LFP'); mats=dir(fullfile(pl,'*_LFP_trialsByRegionAndPhase.mat'));
  pe=fullfile(BASE,subj,'02_Kinematics','Events');
  for f=1:numel(mats)
    M=load(fullfile(pl,mats(f).name)); LP=M.LFP_phases; nt=size(LP,1);
    % matching kinematic file (same block index)
    kf=dir(fullfile(pe,sprintf('*EegPcsEmgKin-%d_kinematic_block.mat',f)));
    if isempty(kf), continue; end
    KB=load(fullfile(kf(1).folder,kf(1).name)); ev=KB.events; vel=double(KB.kinematic_block.velocity); vel=vel(:).';
    % rest beta (block baseline)
    bp=@(x,lo,hi) bandpower(double(x(:))-mean(x),fs,[lo hi]);
    rb=bp(M.baseline_LFP(ch,:),13,30);
    for t=1:min(nt,numel(ev))
      A=ev(t).A;B=ev(t).B;C=ev(t).C;D=ev(t).D;E=ev(t).E;F=ev(t).F;
      if any(isnan([A B C D E F])), continue; end
      reach_dur=(C-A)/fs; grasp_dur=(D-C)/fs; pull_dur=(F-D)/fs; mvt_dur=(F-A)/fs;
      vA=max(1,round(A)); vC=min(numel(vel),round(C)); vD=max(1,round(D)); vF=min(numel(vel),round(F));
      prv=NaN; ppv=NaN;
      if vC>vA, prv=max(vel(vA:vC)); end
      if vF>vD, ppv=max(vel(vD:vF)); end
      % neural per-phase beta (contra STN), dB re block rest
      reach_b=bp(LP{t,2}(ch,:),13,30); grip_b=bp([LP{t,3}(ch,:) LP{t,4}(ch,:)],13,30);
      reach_lb=bp(LP{t,2}(ch,:),13,20); reach_hb=bp(LP{t,2}(ch,:),20,30);
      rlb=bp(M.baseline_LFP(ch,:),13,20); rhb=bp(M.baseline_LFP(ch,:),20,30);
      rows(end+1,:)={subj,f,t,reach_dur,grasp_dur,pull_dur,mvt_dur,prv,ppv,...
        10*log10(reach_b/rb),10*log10(grip_b/rb),10*log10(reach_lb/rlb),10*log10(reach_hb/rhb)}; %#ok<SAGROW>
    end
  end
  fprintf('%s done\n',subj);
end
T=cell2table(rows,'VariableNames',hdr);
writetable(T,fullfile(fileparts(mfilename('fullpath')),'trial_table.csv'));
fprintf('Saved %d trials to trial_table.csv\n',height(T));
