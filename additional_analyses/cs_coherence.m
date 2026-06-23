function [C, f_half] = cs_coherence(S, f)
% CS_COHERENCE  Magnitude-squared coherence from a CSD matrix.
%
%   [C, f_half] = cs_coherence(S, f)
%
%   INPUT
%     S  [m x m x nfft]  two-sided CSD (from cs_csd_multitaper)
%     f  [1 x nfft]      two-sided frequency axis
%
%   OUTPUT
%     C       [m x m x nf]  magnitude-squared coherence in [0,1], one-sided
%                           (0..Nyquist). C(i,j,:) = |S_ij|^2/(S_ii S_jj).
%     f_half  [1 x nf]      one-sided frequency axis (0..Nyquist)

    m    = size(S, 1);
    nfft = size(S, 3);
    nf   = floor(nfft/2) + 1;
    f_half = f(1:nf);

    C = zeros(m, m, nf);
    for i = 1:m
        Sii = real(squeeze(S(i, i, 1:nf))).';
        for j = 1:m
            Sjj = real(squeeze(S(j, j, 1:nf))).';
            Sij = squeeze(S(i, j, 1:nf)).';
            denom = Sii .* Sjj;
            denom(denom <= 0) = eps;
            C(i, j, :) = (abs(Sij).^2) ./ denom;
        end
    end
    C = min(max(C, 0), 1);
end
