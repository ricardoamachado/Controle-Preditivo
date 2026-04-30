function x_saida = diff_n(x_entrada,n)
% Calcula o vetor composto por termos na forma x[i+n] - x[i].
x_saida = x_entrada(1+n:end) - x_entrada(1:end-n);