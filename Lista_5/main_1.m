clear all
close all
clc

l = 2;
m = 200;
E = 0.6e9;
I = 4.17e-5;

M = [m 0 0;
    0 m 0;
    0 0 m];
K = E * I / l^3 * [9/64 1/6 13/192;
            1/6 1/3 1/6;
            13/192 1/6 9/64];


[x,a] = eig(K);

w = sqrt(a/200);

x_c = zeros(4,3);

x_c(2:end,1) = x(:,1) / x(1,1);
x_c(2:end,2) = x(:,2) / x(1,2);
x_c(2:end,3) = x(:,3) / x(1,3);

% Show a figure with a subplot for each column of x_c
figure(1)

% Find overall min and max values for consistent y-limits
y_min = -2;
y_max = 2;

subplot(3,1,1)
plot(x_c(:,1), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
title('{\phi_1}');
grid on;
xlim([1, 4]);
ylim([y_min, y_max]);

subplot(3,1,2)
plot(x_c(:,2), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
title('{\phi_2}');
grid on;
xlim([1, 4]);
ylim([y_min, y_max]);

subplot(3,1,3)
plot(x_c(:,3), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
title('{\phi_3}');
grid on;
xlim([1, 4]);
ylim([y_min, y_max]);

% saving the figure
aveas(gcf, 'modes_1.png');