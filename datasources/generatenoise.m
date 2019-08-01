function [ inp, ts ] = generatenoise( N, freq, length_ms )
%GENERATENOISE - Generate noise with a certain fequency
%   Generate unifrom noise for N neurons over a period of length_ms ms.
%   Assume 1000ms of noise if length not provided. Frequency is in Hz
%   assumes 5 Hz if not provided.
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

inp = [];
ts = [];
diff = num_spikes;

do = true;
while do
    
    inp = [inp, randi([1, N], 1, diff)];
    ts = [ts, randi([1, length_ms], 1, diff)];

    data = unique([ts', inp'], 'rows');
    
    inp = data(:, 2)';
    ts = data(:, 1)';
    diff = num_spikes - numel(inp);
    
    do = diff ~= 0;   % If num_spikes ~= numel(inp)

end

