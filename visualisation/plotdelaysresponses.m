function [ ] = plotdelaysresponses(net, out)
% PLOTDELAYSRESPONSES - Plot the delays and the corresponding response time
% for the output neuron. 


figh = figure('position', [2440, 40, 960, 400]);
subplot(1, 2, 1);
if ~exist('neurons', 'var')
    neurons = [2];
end

if ~exist('asraster', 'var')
    neurons = false;
end

%subplot(2, 1, 1);
if asraster
    plot(out.delays(:, net.N), '.');
    hold on
    for colour = 1 : 23  %Colour code to spike time
        idxs = find(net.pts == colour);
        plot(idxs, out.delays(idxs, net.N), 'o');
    end   
    hold off
else
    plottrace(out.delayst, neurons);
end
title('Pattern synapse delays (sorted)');
axis([0 500 0 20]);

% subplot(2, 1, 2);
% if asraster
%     plot(out.variance(:, net.N), '.');
% else
%     plottrace(out.vart, neurons);
% end
% title('Variance');
% axis([0 2000 0 10]);

end
