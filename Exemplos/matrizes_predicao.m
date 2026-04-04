clc
clear

%% Definindo o sistema em tempo discreto.
A = [0.5 1; 0 0.5];
B = [1;1];
C = [1 1; -1 -1];
D = [0;0];

%% Calculando matrizes de predição.
A_pred = A;
B_pred = B;
num_predicoes = 5;

%Determinação da matriz C de predição.
C_pred = repmat({C},1,num_predicoes);
C_pred = blkdiag(C_pred{:});

for k=2:num_predicoes
    %Determinação da matriz A de predição.
    A_pred = cat(1,A_pred,A ^ k);
    % Primeira coluna da matriz B de predição.
    B_pred = cat(1,B_pred,A ^ (k-1) * B);
end

%Determinação da matriz B de predição.
coluna1_B_pred = B_pred;
for k=1:num_predicoes-1
    %Desloca e gera nova coluna da matriz B.
    nova_col_B_pred = circshift(coluna1_B_pred,k*length(B));
    nova_col_B_pred(1:k*length(B)) = 0;
    B_pred = cat(2,B_pred,nova_col_B_pred);
end

%% Cálculo manual da matriz A de predição para N = 5.
A_pred_esperado = [A;A^2;A^3;A^4;A^5];

assert(isequal(A_pred,A_pred_esperado))

disp(A_pred)

%% Cálculo manual da matriz B de predição para N = 5.
B_pred_esperado = [
    B zeros(size(B)) zeros(size(B)) zeros(size(B)) zeros(size(B))
    A*B B zeros(size(B)) zeros(size(B)) zeros(size(B))
    A^2*B A*B B zeros(size(B)) zeros(size(B));
    A^3*B A^2*B A*B B zeros(size(B));
    A^4*B A^3*B A^2*B A*B B
    ];

assert(isequal(B_pred_esperado, B_pred))
disp(B_pred)
