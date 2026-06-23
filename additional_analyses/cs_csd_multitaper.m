function [S, f] = cs_csd_multitaper(X, fs, NW, K)
% CS_CSD_MULTITAPER  Multitaper cross-spectral density matrix.
%
%   [S, f] = cs_csd_multitaper(X, fs, NW, K)
%
%   Estimates the cross-spectral density (CSD) matrix of a set of signals
%   using Slepian (DPSS) multitaper spectral estimation, averaged across
%   tapers and trials. This single CSD object is the common basis for all
%   downstream connectivity measures (coherence, partial coherence and
%   non-parametric Granger causality) so that every metric is internally
%   consistent.
%
%   INPUT
%     X    [m x n_samp x n_trials]  signals (m nodes). Each (node,trial)
%                                   segment is mean-removed before tapering.
%     fs   scalar                   sampling rate (Hz)
%     NW   scalar                   time-bandwidth product (default 3)
%     K    scalar                   number of tapers (default 2*NW-1)
%
%   OUTPUT
%     S    [m x m x nfft]  TWO-SIDED CSD over bins f = 0 .. fs*(nfft-1)/nfft.
%                          Hermitian per frequency. The full (two-sided)
%                          spectrum is returned because Wilson factorisation
%                          requires it; coherence/partial coherence simply
%                          use the 0..Nyquist half (see cs_*_coherence).
%     f    [1 x nfft]      frequency axis (Hz), two-sided (0..fs).
%
%   The CSD is normalised as a proper power spectral density
%   (divided by fs * sum(taper.^2) and by K*n_trials). The absolute scale
%   does not affect coherence (ratio) or Granger (log-ratio); it only sets
%   the physical units of S and of the noise covariance.
%
%   References: Thomson (1982); Mitra & Pesaran (1999).

    if nargin < 3 || isempty(NW), NW = 3; end
    if nargin < 4 || isempty(K),  K  = 2*NW - 1; end
    K = max(1, round(K));

    [m, n_samp, n_trials] = size(X);
    nfft = n_samp;

    % Slepian tapers: [n_samp x K], each column unit-energy
    tapers = dpss(n_samp, NW, K);            % [n_samp x K]

    S = zeros(m, m, nfft);
    norm_const = fs * n_trials * K;          % PSD normalisation (tapers unit-energy)

    for tr = 1:n_trials
        seg = X(:, :, tr);                   % [m x n_samp]
        seg = seg - mean(seg, 2);            % remove DC per node
        for k = 1:K
            tap = tapers(:, k).';            % [1 x n_samp], unit energy
            F = fft(seg .* tap, nfft, 2);    % [m x nfft] tapered FFT per node
            for a = 1:m
                Fa = F(a, :);
                % accumulate cross-spectra with all nodes b >= a (Hermitian)
                for b = a:m
                    cs = Fa .* conj(F(b, :));      % [1 x nfft]
                    S(a, b, :) = squeeze(S(a, b, :)).' + cs;
                end
            end
        end
    end

    % Fill lower triangle by Hermitian symmetry and normalise
    for a = 1:m
        for b = a:m
            S(a, b, :) = S(a, b, :) / norm_const;
            if b ~= a
                S(b, a, :) = conj(S(a, b, :));
            end
        end
    end

    f = (0:nfft-1) * (fs / nfft);
end
