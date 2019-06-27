function [ output_folder ] = newoutputfolder(basename)
%% NEWOUTPUTFOLDER - return a string of a new enumerated folder 
%   Used to generate new empty folders to dump results into, defaults to
%   the format outxx, where xx is an integer counting up from 00. 
%   TODO - rest of this comment.

count = 0;
if ~exist('basename', 'var')
    basename = 'out';
end

output_folder = [basename, '00'];
while exist(output_folder, 'dir') == 7
    count = count + 1;
    assert(count < 100, 'ERROR: cannot create new results folder, too many folders already exist');
    output_folder = sprintf('%s%02d', basename, count);
end

mkdir(output_folder);

end