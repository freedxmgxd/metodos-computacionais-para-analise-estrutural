clear all
close all
clc

E = 200e9;
A = 0.0127; 

beta=[45 90 135 116.56505 0 135 45 90 90 0 0]*pi/180;
L = [sqrt(2) 3 sqrt(2) sqrt(5) 1 sqrt(2) sqrt(2) 1 1 1 2];

K = zeros(14,14);

% Matrix de equivalencia

ME = [11 12 9 10
    11 12 1 2
    13 14 1 2
    9 10 1 2
    13 14 7 8
    7 8 3 4
    7 8 5 6
    9 10 13 14
    13 14 3 4
    1 2 3 4
    3 4 5 6];

for n=1:11
    k = A*E/L(n);
    c = cos(beta(n));
    s = sin(beta(n));
    ke= k * [c*c c*s -c*c -c*s
        c*s s*s -c*s -s*s
        -c*c -c*s c*c c*s
        -c*s -s*s c*s s*s];
    disp(ke);
    for i=1:4
        for j=1:4
            c=ke(i,j);
            I=ME(n,i);
            J=ME(n,j);
            K(I,J) = K(I,J) + c;
        end
    end
end

disp(K);

% Condições de contorno

% Forças
R = zeros(14,1);
R(6) = -200000;

% D1 = D2 = D11 = 0
D0 = [1 2 11];

K_reduzido = K;
R_reduzido = R;

K_reduzido(D0,:) = [];
K_reduzido(:,D0) = [];
R_reduzido(D0,:) = [];

D_reduzido = inv(K_reduzido)*R_reduzido;

% Reinserindo os valores de D

D = zeros(14,1);
indices_restantes = setdiff(1:14, D0);
D(indices_restantes) = D_reduzido;

tensoes = zeros(11,1);

for i= 1:11
  ui = D(ME(i,1));
  vi = D(ME(i,2));
  uj = D(ME(i,3));
  vj = D(ME(i,4));
  c = cos(beta(i));
  s = sin(beta(i));
  
  deformation = (uj - ui)*c + (vj - vi)*s;
  tensoes(i) = (A*E/L(i))*deformation;
end

disp(tensoes);
