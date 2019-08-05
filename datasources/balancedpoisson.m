function [inp, ts, patt_inp, patt_ts, offsets] = balancedpoisson( Tp, Df, N, Np, Pf, patt_inp, patt_ts, pattfun )
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
%
% Example usage 
%   [ inp, ts, patt_inp, patt_ts, offsets ] = balancedpoisson(50, 10, 2000, 500, 5, [], [])
%
%   NOTE: pattfun not yet supported
%

if exist('pattfun', 'var')
    assert(isempty(pattfun), 'pattfun not yet supported');
end

if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || ...
     isempty(patt_inp) || isempty(patt_ts)
    [ patt_inp, patt_ts ] = generateuniformpattern( Tp, Np );
end

N_inp = N;
inp = [];
ts = [];
max_delay = 20;
slot_size = Tp + max_delay;
num_slots = floor(1000 / slot_size); % TODO dont forget the left over...

probability = Pf / num_slots;
%selected_slots = find(rand(1, num_slots) < probability);
offsets = [];
for r = 1 : 1
    for i = 1 : num_slots
        slot_offset = (i - 1) * slot_size;

        if rand() < probability    % Pattern is presented here

            %% Make 1:Np for time Tp : Tp + max_delay
            non_patt_time = 1000 - Pf * Tp;
            % new frequency necessary to get Df - Pf spikes in whatever time is
            % remaining. e.g. if pattern 5Hz and 50ms long it uses 250ms, thus
            % if Df is 10hz the other fives spikes now need to happen in 750ms.
            % adj_freq is the new frequency necessary.
            adj_freq = (Df - Pf) / (non_patt_time / 1000);
            adj_spikes = adj_freq * Np * (max_delay / 1000); % generate for max_delay ms only.

            % choose adjusted number spikes probabilistically since it doesn't
            % round nicely.
            remainder = adj_spikes - fix(adj_spikes);
            adj_spikes = floor(adj_spikes);
            if rand() > remainder
                adj_spikes = adj_spikes + 1;
            end

            midinp = randi([1, Np], 1, adj_spikes);
            midts = randi([Tp, slot_size], 1, adj_spikes) + slot_offset; 


            %% Make Np + 1 : N_inp from time 1 : Tp
            exp_spikes = ceil(Df * Tp / 1000 * N_inp);
            adj_spikes = exp_spikes - Np;
            topinp = randi([Np + 1, N_inp], 1, adj_spikes);
            topts = randi([1, Tp], 1, adj_spikes) + slot_offset; 


            %% Make Np + 1 : N_inp from time Tp : Tp + max_delay
            all_patt_offset_expected = (N_inp * Df * (Tp / 1000) - Np) * Pf;
            rest_time_expected = Df * (N_inp - Np) - all_patt_offset_expected;
            adj_freq = rest_time_expected / ((1000 - Pf * Tp) / 1000) / (N_inp - Np);
            %exp_spikes = floor(Df * (N_inp) * (slot_size / 1000));
            adj_spikes = adj_freq * (N_inp - Np) * (max_delay / 1000);

            % choose adjusted slot number spike probabilistically
            remainder = adj_spikes - fix(adj_spikes);
            adj_spikes = floor(adj_spikes);
            if rand() > remainder
                adj_spikes = adj_spikes + 1;
            end

            gapinp = randi([Np + 1, N_inp], 1, adj_spikes);
            gapts = randi([Tp, slot_size], 1, adj_spikes) + slot_offset; 


            %% Join all together
            ts = [ts, patt_ts + slot_offset, midts, topts, gapts];
            inp = [inp, patt_inp, midinp, topinp, gapinp];
            offsets = [offsets, slot_offset];

        else                    % Pattern is NOT presented here

            %% Make 1:Np for time 1 : Tp + max_delay
            non_patt_time = 1000 - Pf * Tp;
            % new frequency necessary to get Df - Pf spikes in whatever time is
            % remaining. e.g. if pattern 5Hz and 50ms long it uses 250ms, thus
            % if Df is 10hz the other fives spikes now need to happen in 750ms.
            % adj_freq is the new frequency necessary.
            adj_freq = (Df - Pf) / (non_patt_time / 1000);
            adj_spikes = adj_freq * Np * (slot_size / 1000); % generate for max_delay ms only.

            % choose adjusted number spikes probabilistically since it doesn't
            % round nicely.
            remainder = adj_spikes - fix(adj_spikes);
            adj_spikes = floor(adj_spikes);
            if rand() > remainder
                adj_spikes = adj_spikes + 1;
            end

            midinp = randi([1, Np], 1, adj_spikes);
            midts = randi([1, slot_size], 1, adj_spikes) + slot_offset; 

            %% non pattern gaps
            all_patt_offset_expected = (N_inp * Df * (Tp / 1000) - Np) * Pf;
            rest_time_expected = Df * (N_inp - Np) - all_patt_offset_expected;
            adj_freq = rest_time_expected / ((1000 - Pf * Tp) / 1000) / (N_inp - Np);
            adj_spikes = adj_freq * (N_inp - Np) * (slot_size / 1000);

            % choose adjusted slot number spike probabilistically
            remainder = adj_spikes - fix(adj_spikes);
            adj_spikes = floor(adj_spikes);
            if rand() > remainder
                adj_spikes = adj_spikes + 1;
            end

            rightinp = randi([Np + 1, N_inp], 1, adj_spikes);
            rightts = randi([1, slot_size], 1, adj_spikes) + slot_offset;

            ts = [ts, midts, rightts];
            inp = [inp, midinp, rightinp];
        end


    end
    %% catch the left overs
    slots_end = num_slots * slot_size;
    extra = 1000 - slots_end;
    exp_spikes = (extra / 1000) * Df * N_inp;
    ts = [ts, randi([slots_end, 1000], 1, exp_spikes)];
    inp = [inp, randi([1, N_inp], 1, exp_spikes)];
    
end

end