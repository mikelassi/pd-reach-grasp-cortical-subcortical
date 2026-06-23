function [PC, f_half] = cs_partial_coherence(S, f)
% CS_PARTIAL_COHERENCE  Partial (conditional) magnitude-squared coherence.
%
%   [PC, f_half] = cs_partial_coherence(S, f)
%
%   Partial coherence between nodes i and j removes the linear contribution
%   of ALL other nodes in the system. It is obtained from the inverse CSD
%   (the spectral precision matrix) P = S^{-1}:
%
%       PC_ij = |P_ij|^2 / (P_ii * P_jj)
%
%   A pair that is coherent only because both are driven by a third node
%   will have high ordinary coherence but low partial coherence. This is the
%   key tool for deciding whether cortico-muscular coupling is DIRECT or
%   routed through the STN (and vice versa).
%
%   INPUT
%     S  [m x m x nfft]  two-sided CSD (from cs_csd_multitaper), m >= 2
%     f  [1 x nfft]      two-sided frequency axis
%
%   OUTPUT
%     PC      [m x m x nf]  partial coherence in [0,1], one-sided
%     f_half  [1 x nf]      one-sided frequency axis (0..Nyquist)
%
%   Reference: Dahlhaus (2000); Halliday et al. (1995) partial coherence.

    m    = size(S, 1);
    nfft = size(S, 3);
    nf   = floor(nfft/2) + 1;
    f_half = f(1:nf);

    PC = zeros(m, m, nf);
    ridge = 1e-10;   % tiny diagonal load for numerical stability of inv

    for k = 1:nf
        Sk = S(:, :, k);
        Sk = (Sk + Sk') / 2;                       % enforce Hermitian
        P  = inv(Sk + ridge * trace(Sk)/m * eye(m));
        dP = real(diag(P));
        for i = 1:m
            for j = 1:m
                denom = dP(i) * dP(j);
                if denom <= 0, denom = eps; end
                PC(i, j, k) = (abs(P(i, j))^2) / denom;
            end
        end
    end
    PC = min(max(PC, 0), 1);
end
