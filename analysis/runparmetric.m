% RUNPARMETRIC - Execute the metric calculations on the cluster


filepath = 'mdrp00';
metric = @percentoffsetscorrect;
% metrics = {@percentoffsetscorrect};
% filepaths = {'mdrp00'};
% 
% for metric = metrics
%     for filepath = filepaths
    
        calculatebatchmetric(filepath, metric)
%     end
% end

