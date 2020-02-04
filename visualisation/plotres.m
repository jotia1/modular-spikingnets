%%
cols = size(values, 2);

figh = figure('position', [2440, 40, 480, 480]);
p0 = plot(var_range, values(:, 1:cols), 'x', 'MarkerEdgeColor',[0.7,0.7,0.7]);
hold on
p1 = plot(var_range, mean(values(:, 1:cols), 2), '-.b', 'DisplayName','Mean', 'LineWidth', 2.5);
p2 = plot(var_range, median(values(:, 1:cols), 2), '--r', 'DisplayName','Median', 'LineWidth',2.5);
axis([-Inf Inf 0 1])
ylabel('Accuracy (%)');

legend([p1, p2], 'Location', 'Best')
