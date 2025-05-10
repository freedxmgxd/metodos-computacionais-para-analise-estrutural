clear all
close all
clc

L_total = 3; % m
n_sections = 20;
L = L_total / n_sections; % m

E_beam = 200e9; % Pa
I_beam = 8 * 10^(-4); % m^4
A_beam = 2 * 10^(-2); % m^2
rho_beam = 7800; % kg/m^3

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

D0 = [1, n_sections * 2 + 1];

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

% Plotar os 3 primeiros modos naturais normalizados

figure;
hold on;
colors = {'b', 'r', 'g'}; % Different colors for each mode
legends = {};

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
    
    % Normalize the mode shape to have maximum amplitude of 1
    max_abs_deformation = max(abs(deformation));
    if max_abs_deformation > 0
        deformation = deformation / max_abs_deformation;
    end

    % Plot each mode with a different color
    plot(X, deformation, colors{mode_index}, 'LineWidth', 2);
    legends{end+1} = ['Modo ' num2str(mode_index) ' - ω = ' num2str(modes(mode_index), '%.2f') ' rad/s'];
end

title('Três Primeiros Modos de Vibração Normalizados');
xlabel('Posição (m)');
ylabel('Amplitude Normalizada');
legend(legends);
grid on;
