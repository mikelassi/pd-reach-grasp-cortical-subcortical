%% triad_group_stats.m
% =========================================================================
% Subject-level group inference on the cortex-STN-muscle triad, replacing the
% trial-DOF coherence CI (pseudoreplication) with an exhaustive subject-level
% sign-flip permutation on CONTRASTS. Coherences are atanh(sqrt(.))-stabilised.
%
%   1. Condition modulation of beta coupling (Movement/Pull vs Baseline) for
%      M1contra-STN, STN-muscle, M1contra-muscle.
%   2. ROUTING (3-node triad partial coherence): does conditioning out the
%      third node reduce coupling?
%        - M1-MUS  coh vs pcoh(cond. out STN)  -> is CMC routed via STN?
%        - STN-MUS coh vs pcoh(cond. out M1)   -> is STN-muscle cortically driven?
% =========================================================================
clear; clc;
BASE='F:\Projects\Parkinson_ReachGrasp\Reprocessing';
R=load(fullfile(BASE,'RESULTS_final','CortexSTNMuscle','CortexSTNMuscle_results.mat'));
RES=R.RES; fs=400;
subs=RES.subjects; ncond=numel(RES.cond); nmus=numel(RES.muscles); S=numel(subs);
nl=RES.node_labels;
iM1 = find(strcmpi(nl,'M1_contra')); if isempty(iM1), iM1=find(contains(lower(nl),'m1')&contains(lower(nl),'contra'),1); end
iSTN= find(strcmpi(nl,'STN'));
iMus= find(ismember(nl,RES.muscles));
bidx= find(strcmp({RES.bands.name},'beta (13-30)')); if isempty(bidx), bidx=find(contains({RES.bands.name},'beta'),1); end
fprintf('Nodes: M1_contra=%d STN=%d muscles=%s | beta band idx=%d\n', iM1,iSTN,mat2str(iMus),bidx);

% exhaustive sign matrix
signs=double(dec2bin(0:2^S-1)-'0')*2-1;   % [2^S x S]
atfun=@(x) atanh(sqrt(min(max(x,0),0.999)));
sf=@(d) deal_sf(d,signs);

% ---- 1. Condition modulation of beta coherence (vs Baseline) ----
cBase=find(strcmpi(RES.cond,'Baseline'));
fprintf('\n=== Condition modulation of beta coherence (sign-flip vs Baseline) ===\n');
pairs={'M1-STN',[iM1 iSTN]; 'STN-MUS',[iSTN NaN]; 'M1-MUS',[iM1 NaN]};
for c=1:ncond
  if c==cBase, continue; end
  for pp=1:size(pairs,1)
    nm=pairs{pp,1}; ij=pairs{pp,2};
    if any(isnan(ij))   % muscle pairs: average over muscles
      other=ij(1); vals=nan(S,1); base=nan(S,1);
      for s=1:S
        cc=arrayfun(@(mm) atfun(RES.coh_band(other,iMus(mm),bidx,c,s)),1:nmus);
        bb=arrayfun(@(mm) atfun(RES.coh_band(other,iMus(mm),bidx,cBase,s)),1:nmus);
        vals(s)=nanmean(cc); base(s)=nanmean(bb);
      end
    else
      vals=atfun(squeeze(RES.coh_band(ij(1),ij(2),bidx,c,:)));
      base=atfun(squeeze(RES.coh_band(ij(1),ij(2),bidx,cBase,:)));
    end
    d=vals-base; d=d(~isnan(d)); [t,p]=sf(d);
    fprintf('  %-9s %-8s: dcoh t=%+.2f p=%.4f (n=%d)\n',RES.cond{c},nm,t,p,numel(d));
  end
end

% ---- 2. Routing via triad 3-node partial coherence ----
fprintf('\n=== ROUTING (triad partial coherence, beta, one-sided coh>pcoh) ===\n');
% triad pair columns: 1=M1-STN, 2=STN-MUS, 3=M1-MUS ; pcoh conditions out the 3rd node
ROUT={3,'STN','M1-MUS coupling routed via STN'; 2,'M1','STN-MUS coupling driven via M1'};
for c=1:ncond
  for rr=1:size(ROUT,1)
    col=ROUT{rr,1}; d=nan(S,1);
    for s=1:S
      dm=nan(1,nmus);
      for mm=1:nmus
        C=RES.triad(c).coh{s,mm}; P=RES.triad(c).pcoh{s,mm};
        if isempty(C), continue; end
        nf=size(C,1); fh=(0:nf-1)*(fs/(2*(nf-1))); bb=fh>=13&fh<=30;
        dm(mm)=atfun(mean(C(bb,col)))-atfun(mean(P(bb,col)));
      end
      d(s)=nanmean(dm);
    end
    d=d(~isnan(d)); [t,p]=sf_oneside(d,signs);
    fprintf('  %-9s %-32s: coh-pcoh t=%+.2f p1=%0.4f (%d/%d>0)\n',RES.cond{c},ROUT{rr,3},t,p,sum(d>0),numel(d));
  end
end

% ---- helpers ----
function [t,p]=deal_sf(d,signs)
  S=numel(d); sg=signs(:,1:S); x=d(:).';
  t=mean(x)/(std(x)/sqrt(S)+1e-12);
  perm=mean(sg.*x,2)./(std(sg.*x,0,2)/sqrt(S)+1e-12);
  p=mean(abs(perm)>=abs(t)-1e-9);
end
function [t,p]=sf_oneside(d,signs)
  S=numel(d); sg=signs(:,1:S); x=d(:).';
  t=mean(x)/(std(x)/sqrt(S)+1e-12);
  perm=mean(sg.*x,2)./(std(sg.*x,0,2)/sqrt(S)+1e-12);
  p=mean(perm>=t-1e-9);
end
