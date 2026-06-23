function make_figures(results_dir)
% MAKE_FIGURES  Render all Cortex-STN-Muscle connectivity figures from the
% saved results .mat. Separated from the heavy computation so figures can be
% regenerated/tuned without re-running the analysis.
%
%   make_figures(results_dir)   % folder containing CortexSTNMuscle_results.mat
%   make_figures()              % auto-locate on F: or H:

    if nargin < 1 || isempty(results_dir)
        for d = {'F:\Projects\Parkinson_ReachGrasp\Reprocessing','H:\Parkinson_ReachGrasp\Reprocessing'}
            cand = fullfile(d{1},'RESULTS_final','CortexSTNMuscle');
            if exist(cand,'dir'), results_dir = cand; break; end
        end
    end
    S = load(fullfile(results_dir,'CortexSTNMuscle_results.mat'));
    RES = S.RES;
    bidx = find(strcmp({RES.bands.name},'beta (13-30)'));
    fmin = 2; fmax = 100;

    fig_network(RES, bidx, results_dir);
    fig_triad_undirected(RES, 'coh',  'Coherence',         fmin, fmax, results_dir);
    fig_triad_undirected(RES, 'pcoh', 'Partial coherence', fmin, fmax, results_dir);
    fig_triad_granger(RES, fmin, fmax, results_dir);
    fig_crosssystem_summary(RES, bidx, results_dir);
    fprintf('Figures written to %s\n', results_dir);
end

% -------------------------------------------------------------------------
function fig_network(RES, bidx, outdir)
    labels = RES.node_labels; n = numel(labels);
    cond = RES.cond; meas = {'coh','pcoh','gc'};
    titles = {'Coherence','Partial coherence','Granger (row \leftarrow col)'};
    fig = figure('Visible','off','Position',[40 40 420*numel(cond) 360*3],'Color','w');
    for mi = 1:3
        % shared colour limit across conditions (off-diagonal only)
        vals = [];
        for c = 1:numel(cond)
            M = RES.GA.(meas{mi})(:,:,bidx,c); M(1:n+1:end) = NaN;
            vals = [vals; M(:)]; %#ok<AGROW>
        end
        cmax = prctile(vals(~isnan(vals)), 99); if ~(cmax>0), cmax = eps; end
        for c = 1:numel(cond)
            ax = subplot(3, numel(cond), (mi-1)*numel(cond)+c);
            M = RES.GA.(meas{mi})(:,:,bidx,c); M(1:n+1:end) = NaN;
            imagesc(ax, M, 'AlphaData', ~isnan(M)); axis(ax,'square');
            set(ax,'Color',[0.85 0.85 0.85],'XTick',1:n,'XTickLabel',labels, ...
                'YTick',1:n,'YTickLabel',labels,'FontSize',7,'TickLabelInterpreter','none');
            xtickangle(ax,45); try, clim(ax,[0 cmax]); catch, end
            colormap(ax, hot); cb=colorbar(ax); cb.FontSize=7;
            if mi==1, title(ax,sprintf('%s\n%s (beta)',cond{c},titles{mi}),'FontSize',9);
            else,     title(ax,sprintf('%s (beta)',titles{mi}),'FontSize',9); end
        end
    end
    sgtitle('Cortex-STN-Muscle network — grand-average beta-band connectivity','FontSize',12);
    save_fig(fig, fullfile(outdir,'network_matrices_beta'));
end

% -------------------------------------------------------------------------
function fig_triad_undirected(RES, fld, ttl, fmin, fmax, outdir)
    mus = RES.muscles; nm = numel(mus); cond = RES.cond; nc = numel(cond);
    col = lines(nc); pud = RES.pair_und;
    fig = figure('Visible','off','Position',[40 40 1500 330*nm],'Color','w');
    for m = 1:nm
        for p = 1:3
            ax = subplot(nm,3,(m-1)*3+p); hold(ax,'on'); grid(ax,'on'); ymax = 0;
            hL = [];
            for c = 1:nc
                A = collect(RES.triad(c).(fld), m, p); if isempty(A), continue; end
                f = RES.triad(c).f; fm = f>=fmin & f<=fmax;
                mu = mean(A,2,'omitnan'); se = std(A,0,2,'omitnan')/sqrt(size(A,2));
                fill(ax,[f(fm) fliplr(f(fm))],[(mu(fm)+se(fm))' fliplr((mu(fm)-se(fm))')], ...
                    col(c,:),'EdgeColor','none','FaceAlpha',0.15,'HandleVisibility','off');
                h = plot(ax,f(fm),mu(fm),'Color',col(c,:),'LineWidth',1.8,'DisplayName',cond{c});
                hL(end+1)=h; %#ok<AGROW>
                ymax = max(ymax, max(mu(fm)+se(fm)));
            end
            ymax = max(ymax, 0.02)*1.15;
            ylim(ax,[0 ymax]); xlim(ax,[fmin fmax]);
            % beta band shading to current ylim
            patch(ax,[13 30 30 13],[0 0 ymax ymax],[0.9 0.85 0.1], ...
                'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
            if strcmp(fld,'coh') && ~isnan(RES.ci_line(1))
                yline(ax, RES.ci_line(1), '--', 'Color',[0.4 0.4 0.4], ...
                    'Label','95% CI','FontSize',7,'HandleVisibility','off');
            end
            if m==1, title(ax,pud{p},'FontSize',10); end
            if p==1, ylabel(ax,sprintf('%s\n%s',mus{m},ttl),'FontSize',9); end
            if m==nm, xlabel(ax,'Frequency (Hz)','FontSize',9); end
            if m==1 && p==3 && ~isempty(hL), legend(ax,hL,'Location','northeast','FontSize',7); end
        end
    end
    sgtitle(sprintf('Triad %s — M1_{contra} / STN / Muscle (grand avg \\pm SEM, n=%d)', ...
        ttl, numel(RES.subjects)),'FontSize',12);
    save_fig(fig, fullfile(outdir,sprintf('triad_%s_spectra',fld)));
end

% -------------------------------------------------------------------------
function fig_triad_granger(RES, fmin, fmax, outdir)
    mus = RES.muscles; nm = numel(mus);
    grp = {[1 2],[3 4],[5 6]}; gttl = {'STN \leftrightarrow M1','STN \leftrightarrow MUS','M1 \leftrightarrow MUS'};
    pdir = RES.pair_dir; col2 = [0.10 0.45 0.80; 0.85 0.30 0.10];
    has_thr = isfield(RES.triad,'gc_thr');
    fig = figure('Visible','off','Position',[40 40 1500 330*nm],'Color','w');
    for m = 1:nm
        for g = 1:3
            ax = subplot(nm,3,(m-1)*3+g); hold(ax,'on'); grid(ax,'on'); ymax=0; hL=[];
            dd = grp{g};
            for ci = 1:2
                d = dd(ci);
                A = collect(RES.triad(1).gc, m, d);   % movement
                if isempty(A), continue; end
                f = RES.triad(1).f; fm = f>=fmin & f<=fmax; mu = mean(A,2,'omitnan');
                h = plot(ax,f(fm),mu(fm),'Color',col2(ci,:),'LineWidth',1.8,'DisplayName',pdir{d});
                hL(end+1)=h; ymax=max(ymax,max(mu(fm))); %#ok<AGROW>
                if has_thr
                    T = collect(RES.triad(1).gc_thr, m, d);
                    if ~isempty(T)
                        thr = mean(T,2,'omitnan');
                        plot(ax,f(fm),thr(fm),'--','Color',col2(ci,:),'LineWidth',0.8,'HandleVisibility','off');
                        ymax=max(ymax,max(thr(fm)));
                    end
                end
            end
            ymax = max(ymax,0.01)*1.15; ylim(ax,[0 ymax]); xlim(ax,[fmin fmax]);
            patch(ax,[13 30 30 13],[0 0 ymax ymax],[0.9 0.85 0.1],'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
            if m==1, title(ax,gttl{g},'FontSize',10); end
            if g==1, ylabel(ax,sprintf('%s\nGranger (Movement)',mus{m}),'FontSize',9); end
            if m==nm, xlabel(ax,'Frequency (Hz)','FontSize',9); end
            if ~isempty(hL), legend(ax,hL,'Location','northeast','FontSize',7); end
        end
    end
    ttl = 'Triad directed Granger causality (Movement)';
    if has_thr, ttl = [ttl ' — dashed = surrogate 95% threshold']; end
    sgtitle(ttl,'FontSize',12);
    save_fig(fig, fullfile(outdir,'triad_granger_spectra'));
end

% -------------------------------------------------------------------------
function fig_crosssystem_summary(RES, ~, outdir)
    % Headline: band-averaged coherence for key cross-system links, per condition,
    % shown for the LOW band (corticomuscular common drive) and the BETA band
    % (cortico-subthalamic). Diagonal/self excluded by construction.
    gi = @(s) find(strcmp(RES.node_labels,s));
    bands_to_show = {'low (2-7)','beta (13-30)'};
    links = {'M1_contra','STN','M1-STN'; 'SMA','STN','SMA-STN'; ...
             'PM_contra','STN','PM-STN'; 'STN','IOD','STN-IOD'; ...
             'STN','Triceps','STN-Tri'; 'STN','Deltoid','STN-Del'; ...
             'M1_contra','IOD','M1-IOD'; 'M1_contra','Triceps','M1-Tri'; ...
             'M1_contra','Deltoid','M1-Del'};
    nl = size(links,1); cond = RES.cond; nc = numel(cond);
    fig=figure('Visible','off','Position',[60 60 1200 760],'Color','w');
    for bi = 1:numel(bands_to_show)
        bidx = find(strcmp({RES.bands.name}, bands_to_show{bi}));
        V = zeros(nl,nc);
        for k=1:nl
            i=gi(links{k,1}); j=gi(links{k,2});
            for c=1:nc, V(k,c) = RES.GA.coh(i,j,bidx,c); end
        end
        ax=subplot(2,1,bi); bar(ax, V); grid(ax,'on');
        set(ax,'XTick',1:nl,'XTickLabel',links(:,3),'FontSize',9); xtickangle(ax,25);
        ylabel(ax,sprintf('%s coherence',bands_to_show{bi}),'FontSize',10);
        yline(ax, RES.ci_line(1), '--k', '95% CI','FontSize',9);
        if bi==1, legend(ax, cond, 'Location','northeast','FontSize',9); end
        if bi==1
            title(ax,sprintf(['Cross-system connectivity (n=%d): low-band = corticomuscular ' ...
                'common drive, beta = cortico-subthalamic'], numel(RES.subjects)),'FontSize',11);
        end
    end
    save_fig(fig, fullfile(outdir,'crosssystem_summary'));
end

% -------------------------------------------------------------------------
function A = collect(cellmat, m, col)
    A=[]; if isempty(cellmat), return; end
    for s=1:size(cellmat,1)
        v=cellmat{s,m}; if isempty(v), continue; end
        A=[A, v(:,col)]; %#ok<AGROW>
    end
end
function save_fig(fig, path)
    print(fig, path, '-dpng','-r200'); savefig(fig,[path '.fig']); close(fig);
end
