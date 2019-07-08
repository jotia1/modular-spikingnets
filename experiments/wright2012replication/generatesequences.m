function [ seqs ] = generatesequences(numinputs, numsequences, seqmax)
%% GENERATESEQUENCES - Generate a sequences as per Wright, Wiles (2012)
%   See paper for additonal details, 
%   Returns a matrix with numsequences rows where each row is a sequence,
%   numinputs long in the range [0, seqmax] inclusive. Sequences will be
%   normalised (shifted such that the first firing unit fires at 0)
%
% Assume numinputs == 10, and numsequences == 1 if not supplied
% assume seqmax == 12 if not supplied

if ~exist('numinputs', 'var')
    numinputs = 10;
end
if ~exist('numsequences', 'var')
    numsequences = 2;
end
if ~exist('seqmax', 'var')
    seqmax = 12;
end

% Note +1 and -1 to account for randi not generating zeros. 
seqs = randi(seqmax + 1, numsequences, numinputs) - 1;
seqs = seqs - repmat(min(seqs, [], 2), 1, size(seqs, 2));

end