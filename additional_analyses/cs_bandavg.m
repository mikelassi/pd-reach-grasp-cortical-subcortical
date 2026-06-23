function B = cs_bandavg(M, f, bands)
% CS_BANDAVG  Average a frequency-resolved connectivity array over bands.
%
%   B = cs_bandavg(M, f, bands)
%
%   INPUT
%     M      [.. x nf]   array whose LAST dimension is frequency, OR
%                        [m x m x nf] connectivity cube.
%     f      [1 x nf]    frequency axis matching the last dim of M
%     bands  struct array with field .lim = [lo hi]  (Hz)
%
%   OUTPUT
%     B      same leading dims as M with last dim = numel(bands)

    nb  = numel(bands);
    sz  = size(M);
    nf  = sz(end);
    lead = sz(1:end-1);
    Mr  = reshape(M, [], nf);              % [prod(lead) x nf]
    Br  = zeros(size(Mr,1), nb);
    for b = 1:nb
        m = f >= bands(b).lim(1) & f <= bands(b).lim(2);
        if any(m)
            Br(:, b) = mean(Mr(:, m), 2);
        end
    end
    B = reshape(Br, [lead, nb]);
end
