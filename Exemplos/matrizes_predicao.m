clc
clear

%Definindo o sistema em tempo discreto.
A = [0.5 1; 0 0.5];
B = [1;1];
C = [1 1; -1 -1];
D = [0;0];

%Calculando matrizes de predição.
A_pred = A;
num_predicoes = 5;
for k=2:num_predicoes
    A_pred = cat(1,A_pred,A ^ k);
end
%Cálculo manual da matriz A de predição para N = 5.
A_pred_esperado = [A;A^2;A^3;A^4;A^5];

assert(isequal(A_pred,A_pred_esperado))

disp(A_pred)
disp(A_pred_esperado)