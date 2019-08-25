%%
figure('position', [2440, 40, 640, 480]);
plot(var_range, values(:, 1:cols), 'xk');
hold on
p1 = plot(var_range, mean(values(:, 1:cols), 2), '-.', 'DisplayName','Mean');
p2 = plot(var_range, median(values(:, 1:cols), 2), '--', 'DisplayName','Median');
axis([-Inf Inf 0 1])
ylabel('Accuracy (%)');

legend([p1, p2], 'Location', 'Best')
