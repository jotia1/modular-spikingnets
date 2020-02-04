function [ figh ] = plotall(net, out, plot_lines, do_save)
%% PLOTALL - Plot a series of frequently used graphs in a tabbed figure
%
%   Parameters:
%       net - A network struct 
%       out - A network output struct
%       plot_lines - Whether output spike lines should be drawn
%

% TODO : Hardcoded...
Tp = 50;

figh = figure('name', sprintf('%s : %d_%d : %.2f', net.output_folder, net.count, net.repeat, out.accuracy));
tg = uitabgroup;

%% Create raster plot
rasttab = uitab(tg, 'Title', 'Raster plot');
axes('Parent', rasttab);

if ~exist('do_save','var')
    do_save = false;
end

if ~exist('plot_lines', 'var')
    plot_lines = numel(find(out.spike_time_trace(:, 2) == net.N)) < 1000;
end

rasterspiketimes(out.spike_time_trace, net.group_sizes(1), 1, plot_lines);
highlightoffsets(out.offsets, Tp);
%visualiseaccuracy(net, out);
    
%% Create voltage plot
if numel(net.voltages_to_save) > 0
    volttab = uitab(tg, 'Title', 'voltage plot');
    axes('Parent', volttab);
    plotvoltage(net, out);
end

%% Create variance and delay raster plot
vdrstab = uitab(tg, 'Title', 'var-del raster');
axes('Parent', vdrstab);
plotvardels(net, out, [1], true);

%% Show Iapp plot
if numel(net.iapp_to_save) > 0
    iapptab = uitab(tg, 'Title', 'Iapp');
    axes('Parent', iapptab);
    highlightoffsets(out.offsets, Tp, 0.5);
    hold on
    plot(out.iappt', 'k')
end

%% Show response time plot
vdrstab = uitab(tg, 'Title', 'resp. time');
axes('Parent', vdrstab);
plotresponsetimes(net, out)

if do_save
    savefig(figh, sprintf('%s_%d_%d.fig', net.output_folder, net.count, net.repeat));
end

end