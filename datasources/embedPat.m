function [ inp, ts, patt_inp, patt_ts, patt_locs ] = embedPat( N, patt_inp, patt_ts, sim_time_sec )
%EMBEDPAT return a pattern embedded in noise, optionally takes a pattern
%   If no pattern is given it will create one.
%   Assume make 1 sec of data.

    patt_length = 50; %[ms]
    freq = 10; % spikes per neuron per second [Hz]
    patt_N = ceil(N / 2);
    %patt_N = 300;
    patt_freq = 10; %hz
    if ~exist('sim_time_sec', 'var')
        sim_time_sec = 1;
    end
    sim_time_ms = sim_time_sec * 1000;
    num_presentations = 5;
    total_patt_time = num_presentations * patt_length;
    fired_in_patt = ceil(patt_length / 1000 * freq * patt_N);
    
    if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || numel(patt_inp) == 0 || numel(patt_ts) == 0
       [patt_inp, patt_ts] = genPatt(patt_length, fired_in_patt); 
       disp('Created pattern');
    end

    %% Build non pattern at 10hz
    nonpatt_N = N - patt_N;
    nonpatt_num_spikes = ceil(sim_time_sec * freq * nonpatt_N);
    ts = sort(randi([1, sim_time_ms], 1, nonpatt_num_spikes));
    inp = randi([1, nonpatt_N], 1, nonpatt_num_spikes) + patt_N;
    
    %% Build pattern noise at (10 - 3 =) 7hz
    % num_presentations is the forced frequency
    time_left = sim_time_ms - (num_presentations * patt_length);
    reduced_freq = patt_freq - num_presentations + 1.6;
    %total_nonpatt_time = sim_time_ms - total_patt_time;
    noisepatt_num_spikes = ceil(sim_time_sec * reduced_freq * fired_in_patt);
    noise_ts = sort(randi([1, sim_time_ms], 1, noisepatt_num_spikes));
    noise_inp = randi([1, fired_in_patt], 1, noisepatt_num_spikes);
    %noise_ts = [];
    %noise_inp = [];
    
    %% Trim anything in patterns space
    spacing = ceil(sim_time_ms / num_presentations);
    patt_locs = 0:spacing:900;
    for i = 1:numel(patt_locs)
        start_loc = patt_locs(i);
        idxs = find(noise_ts > start_loc & noise_ts < start_loc + patt_length);
        noise_ts(idxs) = [];
        noise_inp(idxs) = [];
    end
    
    %% GAH HACK
    section_size = spacing - patt_length;
    num_empty_units = patt_N - fired_in_patt;
    empty_spike_num = spacing / 1000 * freq *  num_empty_units;
    empty_ts = [];
    empty_inp = [];
    for i = 1:numel(patt_locs)
        start_loc = patt_locs(i);
        empty_ts = [empty_ts, sort(randi([1, section_size], 1, empty_spike_num)) + start_loc + patt_length];
        empty_inp = [empty_inp, randi([1, num_empty_units], 1, empty_spike_num)+ fired_in_patt];
    end

        
    patt_inp_reps = repmat(patt_inp, [1, num_presentations]);
    patt_ts_reps = repmat(patt_ts, [1, num_presentations]) + sort(repmat(patt_locs, [1, numel(patt_ts)]));
    
    %% Combine all together
    ts = [ts, noise_ts, patt_ts_reps, empty_ts];
    inp = [inp, noise_inp, patt_inp_reps, empty_inp];
    
    % resort according to ts order (maintain ts and inp relative order)
    [ts, idxs] = sort(ts);
    inp = inp(idxs);
   
%       
%     subplot(2, 1, 1)
%     hist(inp, 20)
%     
%     subplot(2, 1, 2)
%     plot(ts, inp, 'k.')
%     hold on
%     
%     % highlight pattern
%     for i = 1:numel(patt_locs)
%         start_loc = patt_locs(i);
%         patch([start_loc start_loc+patt_length start_loc+patt_length start_loc ], [ 0 0 patt_N patt_N ], [0.6 0.4 0.9], 'FaceAlpha',0.5, 'EdgeColor','none')
%     end
%     hold off
%     

end

