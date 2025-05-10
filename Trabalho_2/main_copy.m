clear all
close all
clc

L_total = 200; % in
n_sections = 2;
L = L_total / n_sections; % m

E_beam = 6.58e6; % psi
I_beam = 100; % in^4
A_beam = 2 * 10^(-2); % m^2
rho_beam = 0.1 / A_beam; % kg/m^3

K = zeros((n_sections + 1) * 2);
M = zeros((n_sections + 1) * 2);

for n = 1:n_sections
    % E = 1
    % I = 1
    % A = 1
    % rho = 420
    % L = 1

    K_e = ((E_beam * I_beam) / L^3) * [12, 6 * L, -12, 6 * L;
                                    6 * L, 4 * L^2, -6 * L, 2 * L^2;
                                    -12, -6 * L, 12, -6 * L;
                                    6 * L, 2 * L^2, -6 * L, 4 * L^2];

    M_e = ((rho_beam * A_beam * L) / 420) * [156, 22 * L, 54, -13 * L;
                                        22 * L, 4 * L^2, 13 * L, -3 * L^2;
                                        54, 13 * L, 156, -22 * L;
                                        -13 * L, -3 * L^2, -22 * L, 4 * L^2];

    for i = 1:4

        for j = 1:4
            I = i + (n - 1) * 2;
            J = j + (n - 1) * 2;
            K(I, J) = K(I, J) + K_e(i, j);
            M(I, J) = M(I, J) + M_e(i, j);
        end

    end

end

% Condicoes de contorno
% D1 = D2 = D3 = 0

% Condições de contorno
K_reduzido = K;
M_reduzido = M;

D0 = [1, 2, n_sections * 2 + 1];

K_reduzido(D0, :) = [];
K_reduzido(:, D0) = [];
M_reduzido(D0, :) = [];
M_reduzido(:, D0) = [];

[D_aux, Lambda_aux] = eig(M_reduzido \ K_reduzido);

[lambda, ind] = sort(diag(Lambda_aux));
Lambda = Lambda_aux(ind, ind);
D_reduzido = D_aux(:, ind);

% % Divide cada valor de D pelo primeiro valor de cada coluna
% for i = 1:size(D_reduzido, 2)
%     if D_reduzido(1, i) ~= 0
%         D_reduzido(:, i) = D_reduzido(:, i) / D_reduzido(1, i);
%     else
%         error('Division by zero detected: D(1, %d) is zero.', i);
%     end
% end

% Reinserindo os valores de D
D_modes = zeros((n_sections + 1) * 2, size(D_reduzido, 2));
remaining_indices = setdiff(1:(n_sections + 1) * 2, D0);
D_modes(remaining_indices, :) = D_reduzido;

modes = sqrt(diag(Lambda))
% % Print frequencias naturais
% fprintf('Frequencias naturais:\n');
% for i = 1:length(modes)
%     fprintf('w_%d = %.4f\n', i, modes(i));
% end

% graphics_toolkit("qt");

% Plotar os 3 primeiros modos naturais
figure;

for mode_index = 1:min(3, size(modes))
    % Reset variables for each mode
    deformation = zeros(1, n_sections * 21); % Pre-allocate for better performance
    X = zeros(1, n_sections * 21);
    aux_i = 0;

    for n = 1:n_sections

        for x_ratio = 0:0.05:1% More points for smoother plot
            aux_i = aux_i + 1;
            x = x_ratio * L;

            % Get the modal displacements for this mode
            v1 = D_modes(n * 2 - 1, mode_index);
            theta1 = D_modes(n * 2, mode_index);
            v2 = D_modes(n * 2 + 1, mode_index);
            theta2 = D_modes(n * 2 + 2, mode_index);

            X(aux_i) = x + (n - 1) * L;
            deformation(aux_i) = v1 + theta1 * x + (-(2 * theta1 + theta2) / L - (3 / L^2) * (v1 - v2)) * x^2 + ((theta1 + theta2) / L^2 + (2 / L^3) * (v1 - v2)) * x^3;
        end

    end

    % Plot each mode in a separate subplot
    subplot(3, 1, mode_index);
    plot(X, deformation, 'LineWidth', 2);
    title(['Modo ' num2str(mode_index) ' - ω = ' num2str(modes(mode_index), '%.2f') ' rad/s']);
    xlabel('Posição (in)');
    ylabel('Amplitude');
    grid on;
end

h = 0.003; % s
total_time = 1; % s
n_steps = round(total_time / h);

alpha = 0.00;
beta = 0.00;

C = alpha * M + beta * K;

R = zeros((n_sections + 1) * 2, n_steps);
R_til = zeros((n_sections + 1) * 2, n_steps);

force = 10000; % lb
node = (n_sections + 1);

for n = 1:(0.2 / h)

    if n <= round((0.1 / h) + 1)
        R(node, n) = force;
        % R(node +1, n) = force * L_total / 2;
    else
        R(node, n) = force - (force / (0.1 / h)) * (n - round((0.1 / h) + 1));
        % R(node +1, n) = (force - (force / (0.1 / h)) * (n - round((0.1 / h) + 1))) * L_total / 2;
    end

end

D = zeros((n_sections + 1) * 2, n_steps);
Ddot = zeros((n_sections + 1) * 2, n_steps);
Dddot = zeros((n_sections + 1) * 2, n_steps);

A_1 = ((4 / h^2) * M + (2 / h) * C);
A_2 = ((4 / h) * M + C);
A_3 = M;

Dddot(:, 1) = inv(M) * (R(:, 1) - C * Ddot(:, 1) - K * D(:, 1));
R_til(:, 1) = R(:, 1) + A_1 * D(:, 1) + A_2 * Ddot(:, 1) + A_3 * Dddot(:, 1);

K_til = (4 / h^2) * M + (2 / h) * C + K;

for n = 2:n_steps

    R_til(:, n) = R(:, n) + A_1 * D(:, n - 1) + A_2 * Ddot(:, n - 1) + A_3 * Dddot(:, n - 1);

    D(:, n) = inv(K_til) * R_til(:, n);

    Ddot(:, n) = (2 / h) * D(:, n) - (2 / h) * D(:, n - 1) - Ddot(:, n - 1);
    Dddot(:, n) = (4 / (h^2)) * D(:, n) - (4 / (h^2)) * D(:, n - 1) - (4 / h) * Ddot(:, n - 1) - Dddot(:, n - 1);
end

figure;

% Determine the maximum absolute value from both arrays to set unified limits
max_force_value = force * 2.2;

% Subplot for force
subplot(2, 1, 1);
plot((0:n_steps - 1) * h, R(node, :), 'LineWidth', 2);
title('Força no nó (n_sections + 1)');
xlabel('Tempo (s)');
ylabel('Força (lb)');
% ylim([-max_force_value max_force_value]); % Set unified y-axis limits
grid on;

subplot(2, 1, 2);
plot((0:n_steps - 1) * h, R_til(node, :), 'LineWidth', 2);
title('Força no nó (n_sections + 1) - R_til');
xlabel('Tempo (s)');
ylabel('Força (lb)');
% ylim([-max_force_value max_force_value]); % Set unified y-axis limits
grid on;

% Plotar numa só figura os deslocamentos, velocidades e aceleraçoes no node (n_sections + 1)
figure;
subplot(3, 1, 1);
plot((0:n_steps - 1) * h, D(node, :), 'LineWidth', 2);
title('Deslocamento no nó (n_sections + 1)');
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
% ylim([-2 2]); % Set unified y-axis limits
grid on;
subplot(3, 1, 2);
plot((0:n_steps - 1) * h, Ddot(node, :), 'LineWidth', 2);
title('Velocidade no nó (n_sections + 1)');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
% ylim([-20 20]); % Set unified y-axis limits
grid on;
subplot(3, 1, 3);
plot((0:n_steps - 1) * h, Dddot(node, :), 'LineWidth', 2);
title('Aceleração no nó (n_sections + 1)');
xlabel('Tempo (s)');
ylabel('Aceleração (m/s²)');
% ylim([-200 200]); % Set unified y-axis limits
grid on;
