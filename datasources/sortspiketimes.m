function [inp, ts] = sortspiketimes(inp, ts)
%% SORTSPIKETIMES - Sort the given spike train (inp, ts) by time
    
[ts, idxs] = sort(ts);
inp = inp(idxs);

end