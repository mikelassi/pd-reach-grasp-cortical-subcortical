%% validate_connectivity.m
% =========================================================================
% Validation of the spectral-connectivity estimators (cs_*.m) on synthetic
% systems with KNOWN ground truth. Run this before trusting any real-data
% result. Each test prints PASS/FAIL against an explicit numeric criterion.
%
%   Test 1  Wilson factorisation reconstructs the spectrum  (S ~= H*Z*H')
%   Test 2  Coherence recovers a known analytic value        (mixed latent)
%   Test 3  Partial coherence removes an INDIRECT link        (X->Y->Z chain)
%   Test 4  Granger causality recovers the correct DIRECTION  (X->Y only)
%
% Author: Michael Lassi  (additional_analyses)
% =========================================================================
clear; clc; rng(7);

fs   = 200;            % Hz
T    = 600;            % samples per trial
ntr  = 120;            % trials
NW   = 3;  K = 5;      % multitaper settings
fprintf('=== Connectivity estimator validation ===\n');
fprintf('fs=%d Hz, %d samples x %d trials, multitaper NW=%g K=%d\n\n', fs, T, ntr, NW, K);

allpass = true;

%% ---- TEST 1: Wilson factorisation of an EXACT theoretical spectrum ----
% A known stable MVAR has an analytic CSD that is exactly minimum-phase
% factorisable. Wilson must (a) reconstruct it to ~machine precision and
% (b) recover the true noise covariance Sigma. This isolates algorithm
% correctness from finite-sample estimation noise.
A1 = [ 0.55 0.00;
       0.30 0.55];          % node1->node2 coupling at lag 1
A2 = [-0.80 0.00;
       0.00 -0.80];
Sig  = [1.0 0.0; 0.0 0.7];
nfft = 512;
Sth  = theoretical_csd({A1, A2}, Sig, fs, nfft);   % [2 x 2 x nfft]
[H, Z, Sfac, info] = cs_wilson(Sth, 1e-12, 500);
recon_err = 0; ssum = 0;
for k = 1:nfft
    recon_err = recon_err + norm(Sth(:,:,k) - Sfac(:,:,k), 'fro');
    ssum      = ssum + norm(Sth(:,:,k), 'fro');
end
rel    = recon_err / ssum;
z_err  = norm(Z - Sig, 'fro') / norm(Sig, 'fro');
pass1  = info.converged && rel < 1e-6 && z_err < 1e-3;
fprintf('TEST 1  Wilson (exact spectrum): converged=%d, recon rel.err=%.2e, Sigma rel.err=%.2e  -> %s\n', ...
    info.converged, rel, z_err, tf(pass1));
allpass = allpass && pass1;

%% ---- TEST 2: Coherence recovers a known analytic value ----
% x = s + nx ; y = s + ny ; all independent white, unit variance.
% Theoretical magnitude-squared coherence = vs^2/((vs+vn)(vs+vn)) = 0.25.
vs = 1; vn = 1;
s  = sqrt(vs) * randn(1, T, ntr);
x  = s + sqrt(vn) * randn(1, T, ntr);
y  = s + sqrt(vn) * randn(1, T, ntr);
Xc = [x; y];
[Sc, fc] = cs_csd_multitaper(Xc, fs, NW, K);
[C, fch]  = cs_coherence(Sc, fc);
coh_xy = squeeze(C(1,2,:));
band = fch > 5 & fch < 90;                 % avoid edges
coh_est = mean(coh_xy(band));
coh_true = vs^2 / ((vs+vn)*(vs+vn));        % = 0.25
pass2 = abs(coh_est - coh_true) < 0.03;
fprintf('TEST 2  Coherence: estimated=%.3f, analytic=%.3f  -> %s\n', ...
    coh_est, coh_true, tf(pass2));
allpass = allpass && pass2;

%% ---- TEST 3: Partial coherence removes an indirect link ----
% Chain X -> Y -> Z (no direct X->Z). Ordinary coherence X-Z is nonzero
% (indirect), but partial coherence X-Z|Y must collapse toward zero, while
% the direct links X-Y and Y-Z stay high.
N = T; Z3 = zeros(3, N, ntr);
for tr = 1:ntr
    ex = randn(1,N); ey = randn(1,N); ez = randn(1,N);
    xx = filter(1,[1 -0.6], ex);             % some autocorrelation
    yy = zeros(1,N); zz = zeros(1,N);
    for t = 2:N
        yy(t) = 0.5*yy(t-1) + 0.8*xx(t-1) + 0.5*ey(t);
        zz(t) = 0.5*zz(t-1) + 0.8*yy(t-1) + 0.5*ez(t);
    end
    Z3(:,:,tr) = [xx; yy; zz];
end
[S3, f3] = cs_csd_multitaper(Z3, fs, NW, K);
[C3, fh3]  = cs_coherence(S3, f3);
[PC3, ~]   = cs_partial_coherence(S3, f3);
bb = fh3 > 5 & fh3 < 90;
coh_XZ = mean(squeeze(C3(1,3,bb)));
pc_XZ  = mean(squeeze(PC3(1,3,bb)));
pc_XY  = mean(squeeze(PC3(1,2,bb)));
pc_YZ  = mean(squeeze(PC3(2,3,bb)));
% indirect link must drop substantially under partialisation, while the
% two direct links remain clearly larger than the partialised indirect one.
pass3 = (pc_XZ < 0.5*coh_XZ) && (pc_XY > pc_XZ) && (pc_YZ > pc_XZ);
fprintf('TEST 3  Partial coherence (X->Y->Z chain):\n');
fprintf('        coh(X,Z)=%.3f  -> partial(X,Z|Y)=%.3f   [indirect should drop]\n', coh_XZ, pc_XZ);
fprintf('        partial(X,Y|Z)=%.3f  partial(Y,Z|X)=%.3f  [direct should stay]\n', pc_XY, pc_YZ);
fprintf('        -> %s\n', tf(pass3));
allpass = allpass && pass3;

%% ---- TEST 4: Granger causality recovers the correct direction ----
% Unidirectional AR: x drives y, y does NOT drive x.
A1g = [ 0.55 0.00;
        0.50 0.55];        % A1(2,1)=0.5 : node1(x) -> node2(y)
A2g = [-0.80 0.00;
        0.00 -0.80];
Sg  = [1 0; 0 1];
Xg  = sim_mvar({A1g, A2g}, Sg, T, ntr);
[Sgc, fg] = cs_csd_multitaper(Xg, fs, NW, K);
[G, fgh, gi] = cs_granger(Sgc, fg);
bbg = fgh > 2 & fgh < 90;
gc_x2y = mean(squeeze(G(2,1,bbg)));   % 1->2  (x->y)  expected LARGE
gc_y2x = mean(squeeze(G(1,2,bbg)));   % 2->1  (y->x)  expected ~0
pass4 = (gc_x2y > 0.05) && (gc_x2y > 5*gc_y2x);
fprintf('TEST 4  Granger direction (true: x->y only):\n');
fprintf('        GC(x->y)=%.4f   GC(y->x)=%.4f   ratio=%.1f  -> %s\n', ...
    gc_x2y, gc_y2x, gc_x2y/max(gc_y2x,eps), tf(pass4));
allpass = allpass && pass4;

%% ---- SUMMARY ----
fprintf('\n=========================================\n');
if allpass
    fprintf('ALL TESTS PASSED — estimators are valid.\n');
else
    fprintf('SOME TESTS FAILED — DO NOT trust real-data output yet.\n');
end
fprintf('=========================================\n');

%% ===== local helpers =====
function s = tf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end

function S = theoretical_csd(A, Sig, fs, nfft)
% Analytic two-sided CSD of MVAR(p): z(t)=sum_k A{k} z(t-k)+e, cov(e)=Sig.
%   H(f) = (I - sum_k A{k} exp(-i 2pi f k/fs))^{-1},  S(f) = H Sig H'.
    m = size(A{1},1); p = numel(A);
    S = zeros(m, m, nfft);
    for n = 1:nfft
        fn = (n-1) * fs / nfft;
        Af = eye(m);
        for k = 1:p
            Af = Af - A{k} * exp(-1i*2*pi*fn*k/fs);
        end
        Hf = Af \ eye(m);
        S(:, :, n) = Hf * Sig * Hf';
    end
end

function X = sim_mvar(A, Sig, T, ntr)
% Simulate a stable MVAR(p): z(t) = sum_k A{k} z(t-k) + e(t), cov(e)=Sig.
% A{k}(i,j) is the influence of node j at lag k on node i.
    m = size(A{1},1); p = numel(A);
    L = chol(Sig, 'lower');
    burn = 200;
    X = zeros(m, T, ntr);
    for tr = 1:ntr
        z = zeros(m, T+burn);
        for t = p+1:T+burn
            acc = L * randn(m,1);
            for k = 1:p
                acc = acc + A{k} * z(:, t-k);
            end
            z(:, t) = acc;
        end
        X(:, :, tr) = z(:, burn+1:end);
    end
end
