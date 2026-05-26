clear;

load("sys_disc.mat")
[B_poly, A_poly] = ss2tf(sys_disc.A,sys_disc.B,sys_disc.C,sys_disc.D);

sys_tf = tf(B_poly,A_poly,sys_disc.Ts,'Variable', 'z^-1');
N_predicao = 10;
N_controle = 8;

% Função de transf. na forma (z^-1)*(B(z^-1))/(A(z^-1))
B_poly = B_poly(2:end); % Retira o z^-1.
% Multiplica A por (1-z^-1)
A_til = conv(A_poly, [1, -1]);

%Cada linha j de E representa o polinômio E_j.
%Cada linha j de F representa o polinômio F_j.
[E, F] = Diofantina(A_til,1,N_predicao);

% Calcula o tamanho de cada linha após a convolução.
tamanho_E = size(E, 2);
tamanho_B = length(B_poly);
tamanho_conv = tamanho_E + tamanho_B - 1;
num_linhas = size(E, 1);
M = zeros(num_linhas, tamanho_conv);
%Calcula o produto M_j = E_j * B_j e armazena em M
for idx = 1:num_linhas
    M(idx, :) = conv(E(idx, :), B_poly);
end

% Separa o vetor com coeficientes da resp. ao degrau.
g_vetor = M(N_predicao,1:N_predicao);
g_vetor = g_vetor';

% Monta a matriz da resposta forçada G.
G = g_vetor;
for i=1:N_controle-1
    %Desloca o g_vetor i unidades e preenche início com zeros.
    g_vetor_i = circshift(g_vetor,i);
    g_vetor_i(1:i) = 0;
    G = cat(2,G,g_vetor_i);
end

disp(M);