clc
clear

%% Definindo o sistema em tempo discreto.
load("sys_disc.mat");
A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;

N_estados = size(A,1);
N_entradas = 1;
N_saidas = 1;
% Aumentando o sistema com a ação integral.
Aaug = [A zeros(N_estados,N_saidas);C*A eye(N_saidas)];
Baug = [B;C*B];
Caug = [zeros(N_saidas, N_estados) eye(N_saidas,N_saidas)];
sys_aug = ss(Aaug,Baug,Caug,0,sys_disc.Ts);

%% Matrizes de ponderação.
Q = eye(N_saidas);
R = eye(N_entradas);

%% Calculando matrizes de predição.
%Horizontes.
N_predicao = 10;
N_controle = 10;
%Determinação da matriz C de predição.
C_pred = repmat({Caug},1,N_predicao);
C_pred = blkdiag(C_pred{:});

A_pred = Aaug;
B_pred = Baug;

for k=2:N_predicao
    %Determinação da matriz A de predição.
    A_pred = cat(1,A_pred,Aaug ^ k);
    % Primeira coluna da matriz B de predição.
    B_pred = cat(1,B_pred,Aaug ^ (k-1) * Baug);
end

%Determinação da matriz B de predição.
coluna1_B_pred = B_pred;
for k=1:N_predicao-1
    %Desloca e gera nova coluna da matriz B.
    nova_col_B_pred = circshift(coluna1_B_pred,k*length(Baug));
    nova_col_B_pred(1:k*length(Baug)) = 0;
    B_pred = cat(2,B_pred,nova_col_B_pred);
end
B_pred = B_pred(1:end,1:(N_controle));
%% Determinação da matrizes Q e R de ponderação.
Q_pred = repmat({Q},1,N_predicao);
Q_pred = blkdiag(Q_pred{:});
R_pred = repmat({R},1,N_controle);
R_pred = blkdiag(R_pred{:});

ref_inicial = 2;
W = repmat(ref_inicial,N_predicao,1);
x_inicial = ones(N_estados+1,1);
%% Solução Analítica - Sem Restrições.
delta_U = (B_pred'*C_pred'*Q_pred*C_pred*B_pred + R_pred) ^ (-1) * (B_pred'*C_pred'*Q_pred'*(W - C_pred*A_pred*x_inicial));
