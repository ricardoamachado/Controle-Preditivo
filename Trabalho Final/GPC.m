clear;

load("sys_disc.mat")
A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;
D = sys_disc.D;
T_sample = sys_disc.Ts;

[B_poly, A_poly] = ss2tf(A,B,C,sys_disc.D);
sys_tf = tf(B_poly,A_poly,sys_disc.Ts,'Variable', 'z^-1');
N_predicao = 10;
N_controle = 8;

% Função de transf. na forma (z^-1)*(B(z^-1))/(A(z^-1))
B_poly = B_poly(2:end); % Retira o z^-1.

% Multiplica A por (1-z^-1) para obter o A~.
A_til = conv(A_poly, [1, -1]);

%Cada linha j de E representa o polinômio E_j.
%Cada linha j de F representa o polinômio F_j.
[E, F] = Diofantina(A_til,1,N_predicao);

% Calcula o tamanho de cada linha após a convolução.
tamanho_E = size(E, 2);
tamanho_B = length(B_poly);
tamanho_conv = tamanho_E + tamanho_B - 1;

M = zeros(N_predicao, tamanho_conv);
% Matriz H multiplica os incrementos de controle passados.
% N linhas e n_b colunas, com n_b grau do polinômio B(z^-1).
H = zeros(N_predicao,tamanho_B-1);
%Calcula o produto M_j = E_j * B_j e armazena em M.
for idx = 1:N_predicao
    M(idx, :) = conv(E(idx, :), B_poly);
    H(idx,:) = M(idx, idx+1:idx+tamanho_B-1);
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