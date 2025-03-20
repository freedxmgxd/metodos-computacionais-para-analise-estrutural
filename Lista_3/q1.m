clear all
close all
clc

L = 1; #m
t_simul = 3000; #s

deltax = 0.1; #m
n_knots = L / deltax;
deltat = 750; #s

Ta = 100;
%Ta = 50;
Tb = 100;
%Tb = 50;
alpha = 0.1 * 10^(-4); #m^2/s
% alpha = 1;

s = alpha * deltat / deltax^2;

% s = 1/6;
% deltat = s * deltax^2 / alpha;

steps = int64(t_simul / deltat);

T = zeros(n_knots + 1, steps + 1);
T(1, :) = Ta;
T(n_knots + 1, :) = Tb;

% Condićoes Iniciais
T(1:1) = 50;
T(n_knots + 1, 1) = 50;

% for i = 1:n_knots + 1
%     T(i, 1) = 100 * cos(pi / 2 * ((i - 1) * deltax));
% end

for n = 1:steps

    for i = 2:n_knots
        T_aux = s * T(i - 1, n) + (1 - 2 * s) * T(i, n) + s * T(i + 1, n);

        T(i, n + 1) = T_aux;
    end

    %T(1, n + 1) = T(2, n + 1);
    %T(n_knots + 1, n + 1) = 0;

end

T_exata = zeros(n_knots + 1, 1);

for i = 1:n_knots + 1
    sum_exata = 0;

    M = 100;

    x = (i - 1) * deltax;

    for m = 1:M
        aux_1 = (400 / (((2 * m) -1) * pi));
        aux_2 = sin(((2 * m) - 1) * pi * x)
        aux_3 = exp(-alpha * ((2 * m) - 1)^2 * pi^2 * t_simul)
        sum_exata += aux_1 * aux_2 * aux_3;
    end

    T_exata(i, 1) = 100 - sum_exata;
end

% Show the plot
x = linspace(0, L, n_knots + 1);
t = linspace(0, t_simul, steps + 1);
T_final = T(:, end);
T_final_exata = T_exata(:, 1);

figure(1)
% Plot final numerical solution with scatter points and line
plot(x, T_final, 'r-', 'LineWidth', 2);
hold on
scatter(x, T_final, 50, 'r', 'filled');
% Plot exact solution with dotted line
plot(x, T_final_exata, 'b:', 'LineWidth', 2);

% Add temperature evolution curves with lighter colors and scatter points
for i = 1:floor((steps + 1) / 5):steps% Plot fewer lines to avoid overcrowding
    plot(x, T(:, i), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    scatter(x, T(:, i), 25, [0.7 0.7 0.7], 'filled');
end

xlabel('x [m]');
ylabel('T [°C]');
legend('T numérico', '', 'T exato', 'Evolução temporal', '');
grid on


% Save the figure with high quality
filename = sprintf('graphTx-s%.3f-dx%.2f.png', s, deltax);
set(gcf, 'Position', [100, 100, 800, 600]); % Larger figure size
print(filename, '-dpng', '-r300'); % 300 dpi resolution

% Check if file exists and delete it (will be overwritten anyway by print)
if exist(filename, 'file')
    delete(filename);
    print(filename, '-dpng', '-r300'); % Print again to ensure it's saved
end
