clear all
close all
clc

L_t = 12;
n_sections = 6;

L = L_t / n_sections;
% L = 1;
E = 210 * 10^9;
Ii = 4 * 10^(-4);

k_m = 10 * 10^3;
w = 10 * 10^3;

K = zeros((n_sections + 1) * 2 + 1);

for n = 1:n_sections
    ke = E * Ii / (L^3);
    % ke = 1
    k = ke * [12, 6 * L, -12, 6 * L;
        6 * L, 4 * L * L, -6 * L, 2 * L * L;
        -12, -6 * L, 12, -6 * L;
        6 * L, 2 * L * L, -6 * L, 4 * L * L];

    for i = 1:4

        for j = 1:4
            I = i + (n - 1) * 2;
            J = j + (n - 1) * 2;
            K(I, J) = K(I, J) + k(i, j);
        end

    end

end

k_mola = [k_m, -k_m;
    -k_m, k_m];

k_i = n_sections + 1;
k_f = n_sections * 2 + 3;

K(k_i, k_i) = K(k_i, k_i) + k_mola(1, 1);
K(k_i, k_f) = K(k_i, k_f) + k_mola(1, 2);
K(k_f, k_i) = K(k_f, k_i) + k_mola(2, 1);
K(k_f, k_f) = K(k_f, k_f) + k_mola(2, 2);

% Fazer usando software comercial de elementos finitos para ganhar ponto extra na prova

Q = zeros((n_sections + 1) * 2 + 1, 1);

w_min = w

for n = 1:n_sections / 2
    w_max = w_min;
    w_min = w_max - w / (n_sections / 2)
    [fy_1, m_1, fy_2, m_2] = rampa(L, w_max, w_min);
    Q(n * 2 - 1) = Q(n * 2 - 1) + fy_1;
    Q(n * 2) = Q(n * 2) + m_1;
    Q(n * 2 + 1) = Q(n * 2 + 1) + fy_2;
    Q(n * 2 + 2) = Q(n * 2 + 2) + m_2;
end

w_max = w_min

for n = n_sections / 2 + 1:n_sections
    w_min = w_max
    w_max = w_min + w / (n_sections / 2);
    [fy_1, m_1, fy_2, m_2] = rampa(L, w_max, w_min);
    Q(n * 2 - 1) = Q(n * 2 - 1) + fy_2;
    Q(n * 2) = Q(n * 2) - m_2;
    Q(n * 2 + 1) = Q(n * 2 + 1) + fy_1;
    Q(n * 2 + 2) = Q(n * 2 + 2) - m_1;
end

% Condições de contorno
K_reduzido = K;
Q_reduzido = Q;

D0 = [1, n_sections * 2 + 1, n_sections * 2 + 3];

K_reduzido(D0, :) = [];
K_reduzido(:, D0) = [];
Q_reduzido(D0, :) = [];

D_reduzido = inv(K_reduzido) * Q_reduzido

% Reinserindo os valores de D

D = zeros((n_sections + 1) * 2 + 1, 1);
indices_restantes = setdiff(1:(n_sections + 1) * 2 + 1, D0);
D(indices_restantes) = D_reduzido;

disp(D);
deformation = [];
X = [];
aux_i = 0;

for n = 1:n_sections

    for x = 0:0.5:L
        aux_i = aux_i + 1;

        v1 = D(n * 2 - 1);
        theta1 = D(n * 2);
        v2 = D(n * 2 + 1);
        theta2 = D(n * 2 + 2);

        X(aux_i) = x+(n-1)*L;
        deformation(aux_i) = v1 + theta1 * x + (-(2 * theta1 + theta2) / L -(3 / L^2) * (v1 - v2)) * x^2 + ((theta1+theta2)/L^2 + (2/L^3)*(v1-v2)) * x^3;

    end

end

% Coordenadas dos nós
x = X;
y = zeros(size(x));

% Coordenadas dos nós deformados
fator_escala = 10; % Exagera a deformação
x_deformado = x;
y_deformado = y + fator_escala * deformation;

% Conectividade das barras

barras = zeros(size(x)(2)-1, 2);

for i = 1:size(barras, 1)
    barras(i, 1) = i;
    barras(i, 2) = i + 1;
end

% Plot da treliça na mesma figura
figure;

% Plot da configuração inicial
for i = 1:size(barras, 1)
    no1 = barras(i, 1);
    no2 = barras(i, 2);
    h_inicial = plot([x(no1), x(no2)], [y(no1), y(no2)], 'b-'); % Linhas azuis para a configuração inicial
    hold on;
end

h_inicial_nos = plot(x, y, 'bo', 'MarkerFaceColor', 'b'); % Nós azuis

% Plot da configuração deformada
for i = 1:size(barras, 1)
    no1 = barras(i, 1);
    no2 = barras(i, 2);
    h_deformada = plot([x_deformado(no1), x_deformado(no2)], [y_deformado(no1), y_deformado(no2)], 'r--'); % Linhas vermelhas tracejadas para a configuração deformada
end

h_deformada_nos = plot(x_deformado, y_deformado, 'ro', 'MarkerFaceColor', 'r'); % Nós vermelhos

title('Treliça - Configuração Inicial e Deformada');
axis equal;
legend([h_inicial, h_deformada], 'Inicial', 'Deformada'); % Adiciona legenda correta
