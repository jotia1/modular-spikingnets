function [ inp, ts, patt_inp, patt_ts ] = repeatingrefinedpattern(Tp, Df, N, patt_inp, patt_ts)
%% REPEATINGREFINEDPATTERN - Create 1 second of normalised network input
%   Normalised meaning that creating histagrams (see and of file) will now
%   show any increase in activity during the pattern or any other time. The
%   idea being to make a more challenging task for neural networks to
%   idetify the pattern. Defaults to 100 neurons in pattern (Np)
%
%   Parameters:
%       Tp - Time of pattern (pattern length) [ms]
%       Df - Default frequency of each neuron [Hz]
%       N - Number of neurons in pattern
%       patt_inp/ts - The neurons in the pattern, related to patt_ts
%
% Example usage 
%   [ inp, ts ] = repeatingrefinedpattern(50, 10, 2000, patt_inp, patt_ts)


% Pattern frequence (Pf) is hardcoded for now
% TODO make this more general in the future.
Pf = 2;

%% make pattern
if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || ...
    isempty(patt_inp) || isempty(patt_ts)
    Np = 100;
    patt_inp = 1:Np;
    patt_ts = sort(randi([1, Tp], 1, Np));
else
    Np = numel(patt_inp);
end


% Make rest that isnt pattern
adj_freq = (Df - Pf) / ((1000 - Pf * Tp) / 1000);
adj_spikes = (Df - Pf) * Np;
midinp = randi([1, Np], 1, adj_spikes);
midts = [randi([Tp, 500], 1, adj_spikes / 2), randi([Tp, 500], 1, adj_spikes / 2) + 500]; 

% Make top of pattern
exp_spikes = ceil(Df * Tp / 1000 * N);
adj_spikes = exp_spikes - Np;
topinp = randi([Np + 1, N], 1, adj_spikes * 2);
topts = [randi([1, Tp], 1, adj_spikes), randi([1, Tp], 1, adj_spikes) + 500]; 

% fill gap
exp_spikes = Df * (N - Np);
adj_spikes = exp_spikes - adj_spikes * 2;
gapinp = randi([Np + 1, N], 1, adj_spikes);
gapts = [randi([Tp, 500], 1, adj_spikes / 2), randi([Tp, 500], 1, adj_spikes / 2) + 500]; 


inp = [patt_inp, midinp, patt_inp, topinp, gapinp];
ts = [patt_ts, midts, patt_ts + 500, topts, gapts];

[inp, ts] = sortspiketimes(inp, ts);

%% To verify the data visually we can plot the histograms
if false
    figure;
    
    subplot(2, 1, 1);
    hist(ts);
    title('Ts');
    
    subplot(2, 1, 2);
    hist(inp);
    title('inp');   
end

end