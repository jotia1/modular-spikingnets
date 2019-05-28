function [ res, missing_fields ] = validatenetwork( net )
%VALIDATENETWORK Summary of this function goes here
%   Detailed explanation goes here
res = true;

fields = {'supervising'};
missing_fields = {};

for i = 1:length(fields)
    f = fields{i};
    if ~isfield(net, f)
        res = false;
        missing_fields{end + 1} = f;
    end
end

end

