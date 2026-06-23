function [G, f_half, info] = cs_granger(S, f, tol, maxiter)
% CS_GRANGER  Pairwise non-parametric (spectral) Granger causality.
%
%   [G, f_half, info] = cs_granger(S, f)
%
%   For every ordered pair of nodes, the corresponding 2x2 sub-block of the
%   CSD is factorised with Wilson's algorithm and Geweke's spectral Granger
%   formula is applied. The result is a directed, frequency-resolved measure
%   of how strongly one node's past predicts another's present.
%
%   DIRECTION CONVENTION
%     G(i, j, :)  =  Granger causality  j --> i   (node j drives node i).
%   Diagonal entries are zero.
%
%   INPUT
%     S        [m x m x nfft]  two-sided CSD (from cs_csd_multitaper)
%     f        [1 x nfft]      two-sided frequency axis
%     tol,maxiter              passed to cs_wilson (optional)
%
%   OUTPUT
%     G       [m x m x nf]  directed GC spectra (>=0), one-sided (0..Nyquist)
%     f_half  [1 x nf]      one-sided frequency axis
%     info    struct        .converged [m x m] factorisation convergence
%
%   Pairwise GC is reported (each pair factorised on its own). Spurious
%   coupling driven by a common third node is addressed separately by
%   partial coherence (cs_partial_coherence). Reference: Geweke (1982);
%   Dhamala, Rangarajan & Ding (2008).

    if nargin < 3, tol = []; end
    if nargin < 4, maxiter = []; end

    m    = size(S, 1);
    nfft = size(S, 3);
    nf   = floor(nfft/2) + 1;
    f_half = f(1:nf);

    G = zeros(m, m, nf);
    info.converged = true(m, m);

    for a = 1:m
        for b = a+1:m
            % 2x2 two-sided sub-CSD for nodes (a,b): x=a, y=b
            Sab = S([a b], [a b], :);
            [H, Z, ~, fi] = cs_wilson(Sab, tol, maxiter);
            info.converged(a, b) = fi.converged;
            info.converged(b, a) = fi.converged;

            Sxx = real(squeeze(Sab(1, 1, 1:nf))).';
            Syy = real(squeeze(Sab(2, 2, 1:nf))).';
            Hxy = squeeze(H(1, 2, 1:nf)).';
            Hyx = squeeze(H(2, 1, 1:nf)).';

            Zxx = Z(1, 1);  Zyy = Z(2, 2);  Zxy = Z(1, 2);

            % y -> x : influence of node b on node a
            term_x = (Zyy - Zxy^2 / Zxx) * abs(Hxy).^2;
            gyx = log( Sxx ./ max(Sxx - term_x, eps) );

            % x -> y : influence of node a on node b
            term_y = (Zxx - Zxy^2 / Zyy) * abs(Hyx).^2;
            gxy = log( Syy ./ max(Syy - term_y, eps) );

            gyx = max(real(gyx), 0);
            gxy = max(real(gxy), 0);

            G(a, b, :) = reshape(gyx, 1, 1, nf);   % b -> a
            G(b, a, :) = reshape(gxy, 1, 1, nf);   % a -> b
        end
    end
end
