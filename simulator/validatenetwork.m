function [ isvalid, missing_fields ] = validatenetwork( net )
%VALIDATENETWORK Summary of this function goes here
%   Detailed explanation goes here
isvalid = true;

fields = {'supervising', 'lateral_inhibition'};
missing_fields = {};

for i = 1:length(fields)
    f = fields{i};
    if ~isfield(net, f)
        isvalid = false;
        missing_fields{end + 1} = f;
    end
end

end

