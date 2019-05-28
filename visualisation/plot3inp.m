function [] = plot3inp( out )


%plot(vt(4, :)')
subplot(4, 2, 1)
%plot(debug(:, 1:out.group_sizes(1)));
plot(squeeze(out.delayst(4, :, :))');
title('delays for neuron 4')

subplot(4, 2, 2);
plot(squeeze(out.delayst(5, :, :))');
title('Neuron delays for N5');
%legend({'N1', 'N2', 'N3'});

subplot(4, 2, 3)
plot(squeeze(out.vart(4, :, :))');
title('Variance N4');

subplot(4, 2, 4);
plot(squeeze(out.vart(5, :, :))');
%set(gca, 'Color', 'k')
title('variances of N5');

subplot(2, 2, 3);
%plot(out.vt( net.group_sizes(1) + 1:end, :)')
plot(out.vt(4:5, :)');
hold on
%plot(squeeze(out.v_threst(4, :, :))');
%plot(squeeze(out.v_threst(5, :, :))');
title('output response');
legend({'N4', 'N5', 'n4Thres', 'n5thresh'});
hold off

subplot(2, 2, 4);
plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 4, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 4, 1), 500), '.')
hold on
plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), 500), '.')
hold off
%axis([-Inf Inf 5 14]);
title('output spike time');
legend({'N4', 'N5'});


end