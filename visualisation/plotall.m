function [ ] = plotall(net, out, plot_lines)
%% PLOTALL - Plot a series of frequently used graphs in a tabbed figure
%
%   Parameters:
%       net - A network struct 
%       out - A network output struct
%       plot_lines - Whether output spike lines should be drawn
%

% TODO : Hardcoded...
Tp = 50;

figure;
tg = uitabgroup;

%% Create raster plot
rasttab = uitab(tg, 'Title', 'Raster plot');
axes('Parent', rasttab);

if ~exist('plot_lines', 'var')
    plot_lines = numel(find(out.spike_time_trace(:, 2) == 2001)) < 1000;
end

rasterspiketimes(out.spike_time_trace, 2000, 1, plot_lines);
highlightoffsets(out.offsets, Tp);
    
%% Create voltage plot
if numel(net.voltages_to_save) > 0
    volttab = uitab(tg, 'Title', 'voltage plot');
    axes('Parent', volttab);
    plot(out.vt');
    hold on
    plot( get( gca, 'Xlim' ), [net.v_thres net.v_thres], '--k', 'LineWidth',2)
    highlightoffsets(out.offsets, Tp);
end

%% Create variance and delay raster plot
if numel(net.variance_to_save) > 0 && numel(net.delays_to_save) > 0
    vdrstab = uitab(tg, 'Title', 'var-del raster');
    axes('Parent', vdrstab);
    plotvardels(net, out, [1], true);
end

%% Show Iapp plot
if numel(net.iapp_to_save) > 0
    iapptab = uitab(tg, 'Title', 'Iapp');
    axes('Parent', iapptab);
    highlightoffsets(out.offsets, Tp, 0.5);
    hold on
    plot(out.iappt', 'k')
end


end