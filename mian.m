
n_sections = 6;
K=zeros((n_sections+1)*2+1);

L_t = 12;

L = L_t / n_sections;
L = 1;
E = 210*10^9;
I = 4*10^(-4);

k_m = 10 * 10^3;
w = 1;

for n=1:n_sections
    ke = E*I / L^3;
    k = [12  6*L     -12     6*L
            6*L 4*L*L   -6*L    2*L*L
            -12 -6*L    12      -6*L
            6*L 2*L*L   -6*L    4*L*L];
    for i=1:4
        for j=1:4
            I = i + (n-1)*2;
            J = j + (n-1)*2;
            K(I, J) = K(I, J) + k(i, j);
        end
    end
end

k_mola = [  k_m     -k_m
            -k_m    k_m];

K(7,7) = K(7,7) + k_mola(1,1);
K(7,15) = K(7,15) + k_mola(1,2);
K(15,7) = K(15,7) + k_mola(2,1);
K(15,15) = K(15,15) + k_mola(2,2);



