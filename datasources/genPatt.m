function [ patt_inp, patt_ts ] = genPatt( patt_length, num_spikes )
%GENPATT Summary of this function goes here
%   Detailed explanation goes here

    patt_inp = 1:num_spikes;
    patt_ts = sort(randi([1, patt_length], 1, num_spikes));
    
    %[patt_ts, idxs] = sort(patt_ts);
    %patt_inp = patt_inp(idxs);

end

