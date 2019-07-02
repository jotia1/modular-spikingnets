function [alldata] = plotparamsurface(metric)

folder = 'res';

%% Parameters to sweep
r_fgi = [ 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1 ];
r_a1 =  [ 1, 3, 5, 7 ];
r_a2 = [ 1, 3, 5, 7 ];
r_b1 = [ 1, 3, 5, 7 ];
r_b2 = [ 1, 3, 5, 7 ];

%num_spikes = zeros(4, 1);

alldata = zeros(numel(r_fgi), ...
                numel(r_a1), ...
                numel(r_a2), ...
                numel(r_b1), ...
                numel(r_b2));



%% Run through each variable
for i_b2 = 1 : numel(r_b2)
    for i_b1 = 1 : numel(r_b1) 
        for i_a2 = 1 : numel(r_a2)       
            for i_a1 = 1 : numel(r_a1)
                for i_fgi = 1 : numel(r_fgi)
                    
                    b2 = r_b2(i_b2); b1 = r_b1(i_b1);
                    a2 = r_a2(i_a2); a1 = r_b2(i_a1);
                    fgi = r_fgi(i_fgi);

                    cvt_fgi = sprintf('%.1f', fgi);
                    cvt_fgi(2) = '-';
                    filename = sprintf('%s/%d_%d_%d_%d_%s', folder, a1, a2, b1, b2, cvt_fgi);
                    load(filename, 'net', 'out');
                    
                    value = metric(net, out);
                    
                    alldata(i_fgi, i_a1, i_a2, i_b1, i_b2) = value;
                    
                    %num_spikes(i_fgi) = num_spikes(i_fgi) + totalspikes(net, out);
                    
                end
            end
        end
    end
end

%% Plot results
figure;
num_variables = 4;
steps = 4;

%plottype = @heatmap;
plottype = @(data) imagesc(data);

for i = 1 : steps
    subplot(num_variables, num_variables, 1 + (i-1) * steps);
    plottype(squeeze(alldata(:, :, i, i, i)));
    colorbar;
    title(sprintf('a1-%d', r_a1(i)));

    subplot(num_variables, num_variables, 2 + (i-1) * steps);
    plottype(squeeze(alldata(:, i, :, i, i)));
    colorbar;
    title(sprintf('a2-%d', r_a2(i)));

    subplot(num_variables, num_variables, 3 + (i-1) * steps);
    plottype(squeeze(alldata(:, i, i, :, i)));
    colorbar;
    title(sprintf('b1-%d', r_b1(i)));

    subplot(num_variables, num_variables, 4 + (i-1) * steps);
    plottype(squeeze(alldata(:, i, i, i, :)));
    colorbar;
    title(sprintf('b2-%d', r_b2(i)));
end






% for i = 1 : num_variables
%     for j = 1 : num_variables
%         subplot(num_variables, num_variables, i*j);
%         
%         rows = reshape(repmat(1 : steps, steps, 1), [], 1);
%         cols = repmat(1:steps, step, 1);
%         
%         ind2sub(alldata
%         %plot(alldata(
%         
%     end
% end



end