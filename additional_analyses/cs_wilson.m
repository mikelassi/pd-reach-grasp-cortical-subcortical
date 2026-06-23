function [H, Z, Sfac, info] = cs_wilson(S, tol, maxiter)
% CS_WILSON  Wilson's minimum-phase spectral matrix factorisation.
%
%   [H, Z, Sfac, info] = cs_wilson(S, tol, maxiter)
%
%   Factorises a two-sided CSD matrix S(f) into a minimum-phase transfer
%   function H(f) and a (real, constant) noise covariance Z such that
%
%       S(f) = H(f) * Z * H(f)'      for every frequency f.
%
%   This is the engine for non-parametric (data-driven) Granger causality:
%   it recovers the equivalent linear generative model directly from the
%   spectrum, with NO autoregressive model-order selection.
%
%   INPUT
%     S        [m x m x N]  two-sided CSD over f = 0 .. fs*(N-1)/N (N even),
%                           Hermitian per frequency (from cs_csd_multitaper).
%     tol      scalar       convergence tolerance (default 1e-9)
%     maxiter  scalar       maximum iterations (default 200)
%
%   OUTPUT
%     H     [m x m x N]  minimum-phase transfer function, H(:,:,1) ~ I
%     Z     [m x m]      noise covariance (real, symmetric, positive-def)
%     Sfac  [m x m x N]  reconstructed spectrum H*Z*H' (for validation)
%     info  struct       .converged (logical), .iter, .err
%
%   Reference: Wilson (1972) SIAM J. Appl. Math.; Dhamala, Rangarajan &
%   Ding (2008) NeuroImage 41:354.

    if nargin < 2 || isempty(tol),     tol = 1e-9;  end
    if nargin < 3 || isempty(maxiter), maxiter = 200; end

    m = size(S, 1);
    N = size(S, 3);
    I = eye(m);

    % --- Initialisation: psi(f) = chol(gamma0) for all f ---
    gamma  = ifft(S, [], 3);
    gamma0 = real(gamma(:, :, 1));
    gamma0 = (gamma0 + gamma0') / 2;
    % guarantee positive definiteness for Cholesky
    [R, p] = chol(gamma0);
    if p > 0
        gamma0 = gamma0 + (abs(min(eig(gamma0))) + 1e-12) * eye(m);
        R = chol(gamma0);
    end
    h   = R;                         % gamma0 = h' * h
    psi = repmat(h, [1 1 N]);

    g   = zeros(m, m, N);
    err = Inf;
    converged = false;

    for iter = 1:maxiter
        for f = 1:N
            pf = psi(:, :, f);
            g(:, :, f) = (pf \ S(:, :, f)) / pf' + I;
        end
        gp = plus_operator(g, m, N);

        psi_prev = psi;
        for f = 1:N
            psi(:, :, f) = psi(:, :, f) * gp(:, :, f);
        end

        % convergence: relative change of psi
        num = 0; den = 0;
        for f = 1:N
            num = num + norm(psi(:, :, f) - psi_prev(:, :, f), 'fro');
            den = den + norm(psi(:, :, f), 'fro');
        end
        err = num / max(den, eps);
        if err < tol
            converged = true;
            break;
        end
    end

    % --- Recover H and Z from the converged psi ---
    psi0 = real(ifft(psi, [], 3));
    A0   = psi0(:, :, 1);            % zero-lag term of psi
    Z    = A0 * A0';
    Z    = (Z + Z') / 2;

    H    = zeros(m, m, N);
    Sfac = zeros(m, m, N);
    A0inv = inv(A0);
    for f = 1:N
        H(:, :, f)    = psi(:, :, f) * A0inv;
        Sfac(:, :, f) = psi(:, :, f) * psi(:, :, f)';
    end

    info.converged = converged;
    info.iter = iter;
    info.err  = err;
end

% -------------------------------------------------------------------------
function gp = plus_operator(g, m, N)
% Causal part of g: keep positive lags, half the zero lag, drop negative.
    gam  = real(ifft(g, [], 3));        % lag domain, lag0 at index 1
    gamp = gam;
    gamp(:, :, 1) = 0.5 * gam(:, :, 1); % half the instantaneous term
    gamp(:, :, N/2 + 1 : end) = 0;      % zero the negative lags
    gp = fft(gamp, [], 3);
end
