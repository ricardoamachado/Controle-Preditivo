clear;
clc;

load("sys_disc.mat")
step(sys_disc)
Ts = sys_disc.Ts;
N_amostras = 200;
[g_vetor, t_vetor] = step(sys_disc, 0:Ts:N_amostras*Ts);
% Corta o elemento associado com k = 0.
g_vetor(1) = [];
t_vetor(1) = [];
G = g_vetor;
N_controle = 5;

for i=1:N_controle-1
    %Desloca o g_vetor i unidades e preenche início com zeros.
    g_vetor_i = circshift(g_vetor,i);
    g_vetor_i(1:i) = 0;
    G = cat(2,G,g_vetor_i);
end

Gf = [];
for passo=1:length(g_vetor)
    gf_linha = diff_n(g_vetor,passo)';
    gf_linha = [gf_linha zeros(1,passo)];
    Gf = [Gf ; gf_linha];
end