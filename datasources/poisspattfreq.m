function [ inp, ts, patt_inp, patt_ts, offsets ] = poisspattfreq(Tp, Df, N, Np, Pf, patt_inp, patt_ts, pattfun)
%% POISSPATTFREQ - Return a pattern in distractor noise with poisson freq
%   Return a spike train with a pattern embedded in distractor noise and
%   distributed according to a poisson distribution at the given frequency.
%   Generates 1 second worth of data.
%
%   Note: pattfun should be a function handle that takes a pattern in form
%   (inp, ts) and manipulates it in some way, returning a pattern in form
%   (inp, ts) conforming to Tp, Np.
%
%   Example Usage:
%       [ inp, ts, pattinp, pattts, offsets ] = poisspattfreq(50, 20, 2000, 500, 2, [], [])
%

exp_spikes = poissrnd(Df * N);
inp = randi([1, N], 1, exp_spikes);
ts = randi([1 1000], 1, exp_spikes);

patt_filter = ts > 0 & ts < Tp & inp < Np;

% If no pattern provided, make one
if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || ...
     isempty(patt_inp) || isempty(patt_ts)
    patt_inp = inp(patt_filter);
    patt_ts = ts(patt_filter);
end


num_presentations = min(10, poissrnd(Pf));
for i = 1 : num_presentations
    offset = (i-1) * ceil(1000 / num_presentations);
    
    % Delete old spikes
    replaced_idxs = ts >= offset & ts <= (offset + Tp) & inp <= Np;
    inp(replaced_idxs) = [];
    ts(replaced_idxs) = [];
    
    patt_toinsert = patt_inp;
    ts_toinsert = patt_ts;

    if exist('pattfun', 'var') && isa(pattfun, 'function_handle')
        [patt_toinsert, ts_toinsert] = pattfun(patt_toinsert, ts_toinsert);
    end
    
    %Insert pattern in spots
    inp = [inp, patt_toinsert];
    ts = [ts, ts_toinsert + offset];
    
end

offsets = 0 : ceil(1000 / num_presentations) : 999;
if num_presentations == 0
    offsets = [];
end

end