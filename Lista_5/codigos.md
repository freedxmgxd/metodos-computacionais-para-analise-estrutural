# Codigos Lista 5

## Exercicio 1

```octave
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
```

## Exercicio 2

```octave
clear all
close all
clc

M = [1 0 0;
    0 1 0;
    0 0 2];

C = [0.22 -0.01 0;
    -0.01 0.23 -0.02;
    0 -0.02 0.42];
# C = zeros(3,3);


K = [2 -1 0;
    -1 3 -2;
    0 -2 2];

R = [0;
    0;
    1];

modes = sqrt(eig(M \ K));
% Print frequencias naturais
fprintf('Frequencias naturais:\n');
for i = 1:length(modes)
    fprintf('w_%d = %.4f\n', i, modes(i));
end

T = 2 * pi ./ modes;

h = min(T) / 10;
h = 0.3;
t_final = 70;

n_steps = ceil(t_final / h);

x = zeros(3, n_steps + 1);
x_dot = zeros(3, n_steps + 1);
x_ddot = zeros(3, n_steps + 1);
R_line = zeros(3, n_steps + 1);

x(1, 1) = 0;
x(2, 1) = 0;
x(3, 1) = 0;
x_dot(1, 1) = 0;
x_dot(2, 1) = 0;
x_dot(3, 1) = 0;

t = zeros(1, n_steps + 1);
p = zeros(1, n_steps + 1);

x_ddot(:, 1) = inv(M) * (R - C * x_dot(:, 1) - K * x(:, 1));

A_1 = (4 / (h * h)) * M + (2 / h) * C;
A_2 = (4 / h) * M +C;
A_3 = M;

K_line = (4 / h^2) * M + (2 / h) * C + K;

for i = 2:n_steps + 1
    t(i) = (i - 1) * h;

    R_line(:, i) = R(:) + A_1 * x(:, i - 1) + A_2 * x_dot(:, i - 1) + A_3 * x_ddot(:, i - 1);
    x(:, i) = inv(K_line) * R_line(:, i);
    x_dot(:, i) = 2 / h * x(:, i) - 2 / h * x(:, i - 1) - x_dot(:, i - 1);
    x_ddot(:, i) = (4 / (h * h)) * x(:, i) - (4 / (h * h)) * x(:,i -1) - (4 / h) * x_dot(:, i - 1) - x_ddot(:, i - 1);
end

% Plotting all curves in the same figure
figure(1)
plot(t, x(1, :), 'r-', 'LineWidth', 2);
hold on;
plot(t, x(2, :), 'g-', 'LineWidth', 2);
plot(t, x(3, :), 'b-', 'LineWidth', 2);
hold off;
title('Displacement Responses');
xlabel('Time (s)');
ylabel('Displacement');
legend('x_1', 'x_2', 'x_3');
grid on;

```