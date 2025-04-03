clear all
close all
clc

L = 1; #m
t_simul = 3 * 10^(-1); #s

% deltat = 750; #s

Ta = 100;
%Ta = 50;
Tb = 100;
%Tb = 50;
% alpha = 0.1 * 10^(-4); #m^2/s
alpha = 1;

% s = alpha * deltat / deltax^2;

% Create csv table with s,deltat,deltax,error RMS, razão
filename_csv = 'results.csv';

if exist(filename_csv, 'file')
    delete(filename_csv); % Delete the file if it exists
end

% Open the file for writing
fileID = fopen(filename_csv, 'w');
% Write the header
fprintf(fileID, 's,deltat,deltax,error RMS,razão\n');

S = [0.5, 0.3, 1/6];
deltaX = [0.1, 0.05, 0.025, 0.0125];

for i_s = 1:length(S)
    s = S(i_s);
    RMS_antigo = 0;
    RMS_novo = 0;

    for j_deltaX = 1:length(deltaX)
        deltax = deltaX(j_deltaX);
        n_knots = L / deltax;
        deltat = (s * deltax^2) / alpha;

        steps = int64(t_simul / deltat);

        T = zeros(n_knots + 1, steps + 1);
        T(1, :) = Ta;
        T(n_knots + 1, :) = Tb;

        % %Condićoes Iniciais
        % T(1:1) = 50;
        % T(n_knots + 1, 1) = 50;

        for i = 1:n_knots + 1
            T(i, 1) = 100 * cos(pi / 2 * ((i - 1) * deltax));
        end

        for n = 1:steps

            for i = 2:n_knots
                T_aux = s * T(i - 1, n) + (1 - 2 * s) * T(i, n) + s * T(i + 1, n);

                T(i, n + 1) = T_aux;
            end

            T(1, n + 1) = T(2, n + 1);
            T(n_knots + 1, n + 1) = 0;

        end

        T_exata = zeros(n_knots + 1, steps + 1);
        error_sum = 0;

        for n = 1:steps + 1

            for i = 1:n_knots + 1

                x = double(i - 1) * deltax;
                t = double(n - 1) * deltat;

                T_exata(i, n) = 100 * cos((pi / 2) * x) * exp(-(pi^2/4) * t);

                % Calculate the error
                error_sum += (T(i, n) - T_exata(i, n))^2;
            end

        end

        RMS_antigo = double(RMS_novo);
        RMS_novo = sqrt(double(error_sum) / (double(n_knots + 1) * double(steps + 1)));

        % Write the data
        fprintf(fileID, '%.3f,%.3e,%.2f,%.4f,%.2f\n', s, deltat, deltax, RMS_novo, RMS_antigo / RMS_novo);

        % Show the plot
        x = linspace(0, L, n_knots + 1);
        t = linspace(0, t_simul, steps + 1);
        T_final = T(:, end);
        T_final_exata = T_exata(:, end);

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
        grid on;

        % Save the figure with high quality
        filename = sprintf('graphTx-s%.3f-dx%.2f.png', s, deltax);
        set(gcf, 'Position', [100, 100, 800, 600]); % Larger figure size
        print(filename, '-dpng', '-r300'); % 300 dpi resolution

        % Check if file exists and delete it (will be overwritten anyway by print)
        if exist(filename, 'file')
            delete(filename);
            print(filename, '-dpng', '-r300'); % Print again to ensure it's saved
        end

        close(figure(1)); % Close the figure to avoid displaying it multiple times

    end

end

% Close the file
fclose(fileID);
