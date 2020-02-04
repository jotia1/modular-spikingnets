function [  ] = plotvoltage(net, out)

plot(out.vt');
hold on
%plot( get( gca, 'Xlim' ), [net.v_thres net.v_thres], '--k', 'LineWidth',2)
plot(out.v_threst);
visualiseaccuracy(net, out);
%highlightoffsets(out.offsets, Tp);
%% HACK
hold on;
plot( get( gca, 'Xlim' ), [net.v_thres net.v_thres], '--k', 'LineWidth',1)
plot( get( gca, 'Xlim' ), [-62 -62], '--k', 'LineWidth',1)
%plot(movmean(out.debug, 2000) * 4 - 60, 'k');
bigs = out.debug >= 0;
plot(out.debug / 3 - 62, 'k');

end