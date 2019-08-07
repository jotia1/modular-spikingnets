function [inp, ts, patt_inp, patt_ts, offsets] = balancedpoisson( Tp, Df, N, Np, Pf, patt_inp, patt_ts, pattfun, dropout )
%% BALANCEDPOISSON - Generate carefully balanced poisson input
%   Generate an input stream such that there is no increase in activity
%   either through time or for any given afferent such that detecting a
%   patter is maximally difficult.
%
%   Parameters:
%       Tp - Time of pattern (pattern length) [ms]
%       Df - Default frequency of each neuron [Hz]
%       N - Number of neurons
%       Np - Number of neurons in pattern
%       patt_inp/ts - The neurons in the pattern, related to patt_ts
%       pattfun - the function to apply to the pattern
%       dropout - the percentage to be dropped out
%
% Example usage 
%   [ inp, ts, patt_inp, patt_ts, offsets ] = balancedpoisson(50, 10, 2000, 500, 5, [], [], 0.0)
%

% If there are some dropped by patt fun, compensate by this amount
if ~exist('dropout', 'var')
    dropout = 0.0;
end

if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || ...
     isempty(patt_inp) || isempty(patt_ts)
    [ patt_inp, patt_ts ] = generateuniformpattern( Tp, Np );
end

N_inp = N;
exp_spikes_in_presentation = round(Np * (1 - dropout));

inp = [];
ts = [];

%% Calculate base frequencies 
non_patt_time_sec = (1000 - Pf * Tp) / 1000;
bottom_base_freq = (Df - (Pf * (1 - dropout))) / non_patt_time_sec;

exp_N_inp_spikes_Tp = ceil(Df * Tp / 1000 * N_inp);
exp_N_dist_spikes_Tp = exp_N_inp_spikes_Tp - exp_spikes_in_presentation;
exp_total_N_dist_spikes_Tp = exp_N_dist_spikes_Tp * Pf;

N_dist = N_inp - Np;
exp_N_dist_total_spikes = Df * N_dist;
exp_N_dist_total_spikes_npt = exp_N_dist_total_spikes - exp_total_N_dist_spikes_Tp;
N_dist_npt_freq = exp_N_dist_total_spikes_npt / N_dist / non_patt_time_sec;


%% Generate base spike trains
toppoiss = poissrnd(N_dist_npt_freq * N_dist);
top_inp = randi([Np + 1, N_inp], 1, toppoiss);
top_ts = randi([1, 1000], 1, toppoiss);

bottpoiss = poissrnd(Np * bottom_base_freq);
bott_inp = randi([1, Np], 1, bottpoiss);
bott_ts = randi([1, 1000], 1, bottpoiss);

inp = [ inp, top_inp, bott_inp];
ts = [ts, top_ts, bott_ts];


%% Add offsets
offsets = [];
max_delay = 20;
slot_size = Tp + max_delay;
num_slots = floor(1000 / slot_size);
probability = Pf / num_slots;
selected_slots = find(rand(1, num_slots) < probability);

for i = selected_slots
    offset = (i - 1) * slot_size;
    
    filter = ts > offset & ts < offset + Tp;
    ts(filter) = [];
    inp(filter) = [];

    offsetpoiss = poissrnd(exp_N_dist_spikes_Tp);
    offset_inp = randi([Np + 1, N_inp], 1, offsetpoiss);
    offset_ts = randi([1, Tp], 1, offsetpoiss);

    patt_toinsert = patt_inp;
    ts_toinsert = patt_ts;
    if exist('pattfun', 'var') && isa(pattfun, 'function_handle')
        [patt_toinsert, ts_toinsert] = pattfun(patt_toinsert, ts_toinsert);
    end

    inp = [inp, patt_toinsert, offset_inp];
    ts = [ts, ts_toinsert + offset, offset_ts + offset];
    offsets = [offsets; offset];
end


end