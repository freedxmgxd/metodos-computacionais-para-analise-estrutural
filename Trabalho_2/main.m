clear all
close all
clc

% Carregar pacotes necessários
try
    pkg load video;
    disp('Pacote video carregado com sucesso');
catch
    warning('Falha ao carregar o pacote video. Salvamento de vídeo pode não funcionar.');
end

L_total = 3; % m
n_sections = 20;
L = L_total / n_sections; % m

E_beam = 200e9; % Pa
I_beam = 8 * 10^(-4); % m^4
A_beam = 2 * 10^(-2); % m^2
rho_beam = 7800; % kg/m^3

f_cut = 2000; % Hz

total_time = 0.4; % s

% Condições de contorno
D0 = [1, n_sections * 2 + 1];

K = zeros((n_sections + 1) * 2);
M = zeros((n_sections + 1) * 2);

for n = 1:n_sections

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

K_reduzido = K;
M_reduzido = M;
D_aux_reduzido = zeros((n_sections + 1) * 2, 1);

K_reduzido(D0, :) = [];
K_reduzido(:, D0) = [];
M_reduzido(D0, :) = [];
M_reduzido(:, D0) = [];

[D_aux, Lambda_aux] = eig(M_reduzido \ K_reduzido);

[lambda, ind] = sort(diag(Lambda_aux));
Lambda = Lambda_aux(ind, ind);
D_aux_reduzido = D_aux(:, ind);

% Reinserindo os valores de D
D_modes = zeros((n_sections + 1) * 2, size(D_aux_reduzido, 2));
remaining_indices = setdiff(1:(n_sections + 1) * 2, D0);
D_modes(remaining_indices, :) = D_aux_reduzido;

modes = sqrt(lambda)
frequencys = modes / (2 * pi)

[c index] = min(abs(frequencys - f_cut));

frequency = frequencys(index)

h = 2 * pi / (frequency * 10)
% h = 0.003

critical_damping = 0.05; % Critical damping ratio

A = [1 / modes(1), modes(1);
    1 / modes(2), modes(2); ];

b = 2 * critical_damping * [1;
                        1; ];

x = A \ b;

alpha = x(1)
beta = x(2)
C = alpha * M + beta * K;

C_reduzido = C;
C_reduzido(D0, :) = [];
C_reduzido(:, D0) = [];

n_steps = round(total_time / h);

R = zeros((n_sections + 1) * 2, n_steps);

force_0 = -380000; % N
force_1 = -250000;

L_fraction = 1/4; % Fraction of the total length where the force is applied

element_force = floor(n_sections * L_fraction) + 1;

nodes = [(element_force * 2 - 1), n_sections + 1]

for n = 1:(0.2 / h)

    w = force_0 + ((force_1 - force_0) / (0.2 / h)) * (n - 1);
    x_w = L_total * L_fraction;
    a = w / (L_total - x_w);
    w_b = w;

    for i = 1:n_sections

        if i < element_force
            w_a = 0;
            w_b = 0;
        else
            w_a = w - a * ((i - 1) * L - x_w);
            w_b = w - a * (i * L - x_w);
        end

        R(i * 2 - 1, n) += -w_b * L / 2 - 3 * (w_a - w_b) * L / 20;
        R(i * 2, n) += (-w_b * L^2) / 12 - (w_a - w_b) * L^2/30;
        R(i * 2 + 1, n) += -w_b * L / 2 - 7 * (w_a - w_b) * L / 20;
        R(i * 2 + 2, n) += (w_b * L^2) / 12 + (w_a - w_b) * L^2/20;

    end

end

R_til = zeros(size(R, 1), n_steps);
D = zeros((n_sections + 1) * 2, n_steps);
Ddot = zeros((n_sections + 1) * 2, n_steps);
Dddot = zeros((n_sections + 1) * 2, n_steps);

R(D0, :) = [];
R_til(D0, :) = [];
D(D0, :) = [];
Ddot(D0, :) = [];
Dddot(D0, :) = [];
% % Print frequencias naturais
% fprintf('Frequencias naturais:\n');
% for i = 1:length(modes)
%     fprintf('w_%d = %.4f\n', i, modes(i));
% end

A_1 = ((4 / h^2) * M_reduzido + (2 / h) * C_reduzido);
A_2 = ((4 / h) * M_reduzido + C_reduzido);
A_3 = M_reduzido;

Dddot(:, 1) = inv(M_reduzido) * (R(:, 1) - C_reduzido * Ddot(:, 1) - K_reduzido * D(:, 1));
R_til(:, 1) = R(:, 1) + A_1 * D(:, 1) + A_2 * Ddot(:, 1) + A_3 * Dddot(:, 1);

K_til = (4 / h^2) * M_reduzido + (2 / h) * C_reduzido + K_reduzido;

for n = 2:n_steps

    R_til(:, n) = R(:, n) + A_1 * D(:, n - 1) + A_2 * Ddot(:, n - 1) + A_3 * Dddot(:, n - 1);

    D(:, n) = inv(K_til) * R_til(:, n);

    Ddot(:, n) = (2 / h) * D(:, n) - (2 / h) * D(:, n - 1) - Ddot(:, n - 1);
    Dddot(:, n) = (4 / (h^2)) * D(:, n) - (4 / (h^2)) * D(:, n - 1) - (4 / h) * Ddot(:, n - 1) - Dddot(:, n - 1);
end

R_expanded = zeros((n_sections + 1) * 2, n_steps);
R_til_expanded = zeros((n_sections + 1) * 2, n_steps);
D_expanded = zeros((n_sections + 1) * 2, n_steps);
Ddot_expanded = zeros((n_sections + 1) * 2, n_steps);
Dddot_expanded = zeros((n_sections + 1) * 2, n_steps);

R_expanded(remaining_indices, :) = R;
R_til_expanded(remaining_indices, :) = R_til;
D_expanded(remaining_indices, :) = D;
Ddot_expanded(remaining_indices, :) = Ddot;
Dddot_expanded(remaining_indices, :) = Dddot;

% Plotar os 3 primeiros modos naturais normalizados
% Create directory if it doesn't exist
output_dir = 'output_figures';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

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
    legends{end + 1} = ['Modo ' num2str(mode_index) ' - ω = ' num2str(modes(mode_index), '%.2f') ' rad/s'];
end

title('Três Primeiros Modos de Vibração Normalizados');
xlabel('Posição (m)');
ylabel('Amplitude Normalizada');
legend(legends);
grid on;

% Define all figure parameters in a single array
fig_params = [0, 0, 1280, 720];

% Save the figure in the created directory with specified resolution
filename = fullfile(output_dir, 'modos_vibracao');
set(gcf, 'Position', fig_params); % Set figure size
print(filename, '-dpng', '-r100'); % PNG format with specified DPI
disp(['Figure saved to: ', filename]);

% Plotar uma coluna de D em uma figura

figure;

for i = 1:n_steps
    plot((0:L_total / n_sections:L_total), transpose(D_expanded(1:2:end, i)), 'k--o', 'LineWidth', 0.1);
    hold on;
end

plot((0:L_total / n_sections:L_total), transpose(D_expanded(1:2:end, 1)), 'r--.', 'LineWidth', 0.5);
title('Envelope do deslocamento');
xlabel('Posição (m)');
ylabel('Deslocamento (m)');
grid on;

% Save the figure in the created directory with specified resolution
filename = fullfile(output_dir, 'envelope_deslocamento');
set(gcf, 'Position', fig_params); % Set figure size
print(filename, '-dpng', '-r100'); % PNG format with specified DPI
disp(['Figure saved to: ', filename]);

figure;
% Subplot for force
subplot(2, 1, 1);
plot((0:n_steps - 1) * h, R_expanded(nodes(1), :), 'k--o', 'LineWidth', 2);
title(['Força no nó ', num2str((nodes(1) + 1) / 2)]);
xlabel('Tempo (s)');
ylabel('Força (N)');
% ylim([-max_force_value max_force_value]); % Set unified y-axis limits
grid on;

subplot(2, 1, 2);
plot((0:L_total / n_sections:L_total), transpose(R_expanded(1:2:end, 1)), 'ko', 'LineWidth', 0.5);
title(['Força no t=0']);
xlabel('Posição (m)');
ylabel('Força (N)');
% ylim([-max_force_value max_force_value]); % Set unified y-axis limits
grid on;

% Save the figure in the created directory with specified resolution
filename = fullfile(output_dir, 'forca');
set(gcf, 'Position', fig_params); % Set figure size
print(filename, '-dpng', '-r100'); % PNG format with specified DPI
disp(['Figure saved to: ', filename]);

for node = nodes
    % Plotar numa só figura os deslocamentos, velocidades e aceleraçoes no node (n_sections + 1)
    figure;
    subplot(3, 2, 1);
    plot((0:n_steps - 1) * h, D_expanded(node, :), 'LineWidth', 2);
    title(['Deslocamento no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Deslocamento (m)');
    % ylim([-2 2]); % Set unified y-axis limits
    grid on;
    subplot(3, 2, 3);
    plot((0:n_steps - 1) * h, Ddot_expanded(node, :), 'LineWidth', 2);
    title(['Velocidade no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Velocidade (m/s)');
    % ylim([-20 20]); % Set unified y-axis limits
    grid on;
    subplot(3, 2, 5);
    plot((0:n_steps - 1) * h, Dddot_expanded(node, :), 'LineWidth', 2);
    title(['Aceleração no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Aceleração (m/s²)');
    % ylim([-200 200]); % Set unified y-axis limits
    grid on;
    subplot(3, 2, 2);
    plot((0:n_steps - 1) * h, D_expanded(node + 1, :), 'LineWidth', 2);
    title(['Rotação no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Rotação (rad)');
    % ylim([-2 2]); % Set unified y-axis limits
    grid on;
    subplot(3, 2, 4);
    plot((0:n_steps - 1) * h, Ddot_expanded(node + 1, :), 'LineWidth', 2);
    title(['Velocidade angular no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Velocidade (rad/s)');
    % ylim([-20 20]); % Set unified y-axis limits
    grid on;
    subplot(3, 2, 6);
    plot((0:n_steps - 1) * h, Dddot_expanded(node + 1, :), 'LineWidth', 2);
    title(['Aceleração angular no nó ', num2str((node + 1) / 2)]);
    xlabel('Tempo (s)');
    ylabel('Aceleração (rad/s²)');
    % ylim([-200 200]); % Set unified y-axis limits
    grid on;

    % Save the figure in the created directory
    name = ['deslocamento_forca_' num2str(node)];
    filename = fullfile(output_dir, '', name);
    set(gcf, 'Position', fig_params); % Set figure size
    print(filename, '-dpng', '-r100'); % PNG format with specified DPI
    disp(['Figure saved to: ', filename]);

end

% Create figure for displacement animation
figure;
set(gcf, 'Position', fig_params); % Set figure size for better quality

% Create the main plot axes
plot_handle = plot((0:L_total / n_sections:L_total), transpose(D_expanded(1:2:end, 1)), 'k--o', 'LineWidth', 1.5);
title('Deslocamento ao longo do tempo');
xlabel('Posição (m)');
ylabel('Deslocamento (m)');
grid on;

% Find the maximum displacement for consistent y-axis limits
max_disp = max(max(abs(D_expanded(1:2:end, :))));
ylim([-max_disp * 1.1 max_disp * 1.1]);

% Text display for time
time_text = text(0.05, 0.95, 'Tempo: 0.000 s', 'Units', 'normalized');

% Setup for saving animation
video_filename = fullfile(output_dir, 'displacement_animation.avi');
% Use default AVI format which is guaranteed to work with Octave
video_writer = VideoWriter(video_filename);
disp('Using AVI format for video output');
video_writer.FrameRate = 30; % Frame rate of the output video
open(video_writer);

% Animation parameters
% frame_skip = max(1, floor(n_steps / 100)); % Skip frames to speed up animation if too many steps
pause_time = 0.01; % Time between frames in seconds

% Run the animation
for i = 1:n_steps
    % Update the plot data
    set(plot_handle, 'YData', transpose(D_expanded(1:2:end, i)));

    % Update time display
    set(time_text, 'String', ['Tempo: ' num2str((i - 1) * h, '%.3f') ' s']);

    % Refresh the display
    drawnow;

    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(video_writer, frame);

    % Pause to control animation speed (only when displaying)
    pause(pause_time);
end

% Close the video file
close(video_writer);
disp(['Animation saved to: ', video_filename]);

% Create figure for displacement and force animation
figure;
set(gcf, 'Position', fig_params); % Set figure size for better quality

% Create subplot for displacement
subplot(2, 1, 1);
disp_handle = plot((0:L_total / n_sections:L_total), transpose(D_expanded(1:2:end, 1)), 'k--o', 'LineWidth', 1.5);
title('Deslocamento ao longo do tempo');
xlabel('Posição (m)');
ylabel('Deslocamento (m)');
grid on;

% Find the maximum displacement for consistent y-axis limits
max_disp = max(max(abs(D_expanded(1:2:end, :))));
ylim([-max_disp * 1.1 max_disp * 1.1]);

% Create subplot for force using stem plot instead of plot
subplot(2, 1, 2);
force_handle = stem((0:L_total / n_sections:L_total), transpose(R_expanded(1:2:end, 1)), 'r', 'LineWidth', 1.5, 'MarkerSize', 4);
title('Força nodal equivalente ao longo do tempo');
xlabel('Posição (m)');
ylabel('Força (N)');
grid on;

% Find the maximum force for consistent y-axis limits
max_force = max(max(abs(R_expanded(1:2:end, :))));
ylim([-max_force * 1.1 max_force * 1.1]);

% Text display for time (place it on top subplot)
subplot(2, 1, 1);
time_text = text(0.05, 0.95, 'Tempo: 0.000 s', 'Units', 'normalized');

% Setup for saving animation
video_filename = fullfile(output_dir, 'displacement_force_animation.avi');
% Use default AVI format which is guaranteed to work with Octave
video_writer = VideoWriter(video_filename);
disp('Using AVI format for video output');
video_writer.FrameRate = 30; % Frame rate of the output video
open(video_writer);

% Animation parameters
pause_time = 0.01; % Time between frames in seconds

% Run the animation
for i = 1:n_steps
    % Update displacement plot
    subplot(2, 1, 1);
    set(disp_handle, 'YData', transpose(D_expanded(1:2:end, i)));

    % Update time display
    set(time_text, 'String', ['Tempo: ' num2str((i - 1) * h, '%.3f') ' s']);

    % Update force plot - need special handling for stem
    subplot(2, 1, 2);
    delete(force_handle);
    force_handle = stem((0:L_total / n_sections:L_total), transpose(R_expanded(1:2:end, i)), 'r', 'LineWidth', 1.5, 'MarkerSize', 4);
    ylim([-max_force * 1.1 max_force * 1.1]);
    title('Força nodal equivalente ao longo do tempo');
    xlabel('Posição (m)');
    ylabel('Força (N)');
    grid on;
    % Refresh the display
    drawnow;

    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(video_writer, frame);

    % Pause to control animation speed (only when displaying)
    pause(pause_time);
end

% Close the video file
close(video_writer);
disp(['Animation saved to: ', video_filename]);
