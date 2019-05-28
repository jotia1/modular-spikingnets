function [ inp, ts ] = generatenoise( N, freq, length_ms )
%GENERATENOISE - Generate noise with a certain fequency
%   Generate unifrom noise for N neurons over a period of length_ms ms.
%   Assume 1000ms of noise if length not provided. Frequency is in Hz
%   assumes 5 Hz if not provided.
%
%   NOTE:
%   Duplicates are removed at the end meaning if there are lots of clashes
%   (caused by high frequency in a small period) the statistics may be 
%   quite off. 
%
%   TODO
%       - Replace removed items
%
%   Example usage:
%       [ inp, ts ] = generatenoise( 10, 2, 1000 )

if ~exist('freq', 'var') || isempty(freq)
    freq = 5;
end

if ~exist('length_ms', 'var') || isempty(length_ms)
    length_ms = 1000;
end


length_sec = length_ms / 1000;
num_spikes = ceil(freq * length_sec * N);

inp = randi([1, N], num_spikes, 1);
ts = randi([1, length_ms], num_spikes, 1);

data = unique([ts, inp], 'rows');
inp = data(:, 2);
ts = data(:, 1);

end

