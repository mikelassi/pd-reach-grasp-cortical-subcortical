% check_results.m -- group-level summary of the Cortex-STN-Muscle results
base=''; for d={'F:\Projects\Parkinson_ReachGrasp\Reprocessing','H:\Parkinson_ReachGrasp\Reprocessing'}
    if exist(d{1},'dir'), base=d{1}; break; end, end
R=load(fullfile(base,'RESULTS_final','CortexSTNMuscle','CortexSTNMuscle_results.mat')); RES=R.RES;
fprintf('Subjects (n=%d): %s\n', numel(RES.subjects), strjoin(RES.subjects,', '));
fprintf('Bands: %s\n', strjoin({RES.bands.name},' | '));
fprintf('CI(95%%) coherence line per cond: %s\n', mat2str(RES.ci_line,3));
gi=@(s) find(strcmp(RES.node_labels,s));
lowi=find(strcmp({RES.bands.name},'low (2-7)'));
beti=find(strcmp({RES.bands.name},'beta (13-30)'));
condn={'Mvt','Pull','Base'};

links={'M1_contra','STN';'SMA','STN';'PM_contra','STN'; ...
       'M1_contra','IOD';'M1_contra','Triceps';'M1_contra','Deltoid'; ...
       'STN','IOD';'STN','Triceps';'STN','Deltoid'};
for bb=[lowi beti]
    fprintf('\n=== Grand-avg %s coherence (Mvt/Pull/Base), CI=%.3f  (*=>CI) ===\n', RES.bands(bb).name, RES.ci_line(1));
    for k=1:size(links,1)
        i=gi(links{k,1}); j=gi(links{k,2});
        v=[RES.GA.coh(i,j,bb,1) RES.GA.coh(i,j,bb,2) RES.GA.coh(i,j,bb,3)];
        st=''; if max(v)>RES.ci_line(1), st=' *'; end
        fprintf('  %-10s - %-9s : %.3f / %.3f / %.3f%s\n', links{k,1},links{k,2}, v, st);
    end
end

fprintf('\n=== Triad coherence PEAK over subjects (Movement) low/beta/hibeta ===\n');
for m=1:numel(RES.muscles)
    for p=1:3
        A=gather(RES.triad(1).coh,m,p); f=RES.triad(1).f;
        if isempty(A), continue; end
        mu=mean(A,2,'omitnan');
        bl=f>=2&f<=7; bb2=f>=13&f<=30; bh=f>=30&f<=45;
        fprintf('  %-8s %-8s: low=%.3f@%.0f  beta=%.3f@%.0f  hibeta=%.3f@%.0f\n', ...
            RES.muscles{m}, RES.pair_und{p}, max(mu(bl)),pk(mu,f,bl), max(mu(bb2)),pk(mu,f,bb2), max(mu(bh)),pk(mu,f,bh));
    end
end

% Directed Granger vs surrogate threshold (Movement), low & beta bands
fprintf('\n=== Triad directed Granger (Movement): obs vs surrogate-95 (low | beta) ===\n');
hasT = isfield(RES.triad,'gc_thr');
for m=1:numel(RES.muscles)
    fprintf('-- %s --\n', RES.muscles{m});
    for d=1:6
        A=gather(RES.triad(1).gc,m,d); f=RES.triad(1).f;
        if isempty(A), continue; end
        mu=mean(A,2,'omitnan');
        bl=f>=2&f<=7; bb2=f>=13&f<=30;
        sl=''; sb='';
        if hasT
            T=gather(RES.triad(1).gc_thr,m,d);
            if ~isempty(T), thr=mean(T,2,'omitnan');
                if mean(mu(bl))>mean(thr(bl)), sl='*'; end
                if mean(mu(bb2))>mean(thr(bb2)), sb='*'; end
            end
        end
        fprintf('   %-10s: low=%.4f%s  beta=%.4f%s\n', RES.pair_dir{d}, mean(mu(bl)),sl, mean(mu(bb2)),sb);
    end
end

function A=gather(cm,m,col)
    A=[]; if isempty(cm), return; end
    for s=1:size(cm,1), v=cm{s,m}; if ~isempty(v), A=[A,v(:,col)]; end, end
end
function fp=pk(mu,f,b), fb=f(b); mub=mu(b); [~,i]=max(mub); fp=fb(i); end
