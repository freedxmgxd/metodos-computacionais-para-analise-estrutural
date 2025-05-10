clear all
close all
clc

m = 0.45594;
k = 18;
c = 2*0.2865;

t_final = 6;
n_steps = 100;
dt = t_final / n_steps;

t = zeros(1, n_steps + 1);
p = zeros(1, n_steps + 1);

t_pulse = 1/2 * (2 * pi / 5.236);

u = zeros(1, n_steps + 1);
u_dot = zeros(1, n_steps + 1);
u_ddot = zeros(1, n_steps + 1);
p_line = zeros(1, n_steps + 1);
u(1) = 0;
u_dot(1) = 0;

t(1) = 0;
p(1) = 50 * sin(5.236 * t(1));
u_ddot(1) = (p(1) - k * u(1) - c * u_dot(1)) / m;

h = dt;
a_1 = (4 * m / (h * h)) + (2 * c / h);
a_2 = (4 * m / h) +c;
a_3 = m;
k_line = a_1 + k;

for i = 2:n_steps + 1
    t(i) = (i - 1) * dt;

    if t(i) > t_pulse
        p(i) = 0;
    else
        p(i) = 50 * sin(5.236 * t(i));
    end

    p_line(i) = p(i) + a_1 * u(i - 1) + a_2 * u_dot(i - 1) + a_3 * u_ddot(i - 1);
    u(i) = p_line(i) / k_line;
    u_dot(i) = ((- 2) / h) * u(i - 1) - u_dot(i - 1) + (2 / h) * u(i);
    u_ddot(i) = ((-4) / (h * h)) * u(i - 1) - (4 / h) * u_dot(i - 1)- u_ddot(i - 1) + (4 / (h * h)) * u(i);
end

figure(1)
plot(t, p, 'r', 'LineWidth', 2);
xlabel('t [s]');
ylabel('p(t) [N]');
title('Forcing Function p(t)');
grid on;
 
figure(2)
plot(t, u, 'b', 'LineWidth', 2);
xlabel('t [s]');
ylabel('u(t) [m]');
title('Displacement u(t)');
grid on;

