% Define the vectors
w1 = [0, 1];
w0 = [1, 0];
w_tau = [0.7, 0.7]; % Adjust as needed to show the Pareto optimal vector

% Create the plot
figure;
hold on;
axis equal;
grid on;

% Plot the vectors
quiver(0, 0, w1(1), w1(2), 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5); % w1 in blue
quiver(0, 0, w0(1), w0(2), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % w0 in orange
quiver(0, 0, w_tau(1), w_tau(2), 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % w_tau in green

% Plot the dashed Pareto optimal curve (approximate with an arc)
theta = linspace(0, pi/2, 100);
x_arc = cos(theta);
y_arc = sin(theta);
plot(x_arc, y_arc, 'k--');

% Add text labels
text(w1(1), w1(2), '\bf C-C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(w0(1), w0(2), '\bf S-C', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(w_tau(1), w_tau(2), '\bf Pareto Optimal', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(0, 1, 'hc/||hc||', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(1, 0, 'hs/||hs||', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Set axis limits
xlim([-0.2, 1.2]);
ylim([-0.2, 1.2]);

% Axis labels
xlabel('Real Part');
ylabel('Imaginary Part');

% Title
title('Vector Plot with Pareto Optimal Curve');

hold off;
