clc
clear

%% Definindo o sistema em tempo discreto.
load("sys_disc.mat");
A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;

N_estados = size(A,1);
N_entradas = size(B,2);
N_saidas = size(C,1);
% Aumentando o sistema com a ação integral.
Aaug = [A zeros(N_estados,N_saidas);C*A eye(N_saidas)];
Baug = [B;C*B];
Caug = [zeros(N_saidas, N_estados) eye(N_saidas,N_saidas)];
sys_aug = ss(Aaug,Baug,Caug,0,sys_disc.Ts);

%% Matrizes de ponderação.
Q = eye(num_saidas);
R = eye(num_entradas);

%% Calculando matrizes de predição.
A_pred = Aaug;
B_pred = Baug;
N_predicao = 10;
%TODO: Continuar montagem da predição.