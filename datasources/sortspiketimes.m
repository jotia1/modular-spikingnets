function [inp, ts] = sortspiketimes(inp, ts)
    
[ts, idxs] = sort(ts);
inp = inp(idxs);

end