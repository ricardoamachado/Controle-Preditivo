clear
clc

%Definindo o sistema em tempo discreto.
A = [0.5 1; 0 0.5];
B = [1;1];
C = [1 1; -1 -1];
D = [0;0];
sys = ss(A,B,C,D,-1);

num_estados = size(A,1);
num_entradas = size(B,2);
num_saidas = size(C,1);

%Implementação da ação integral.
Aaug = [A zeros(num_estados,num_saidas);C*A eye(num_saidas)];
Baug = [B;C*B];
Caug = [zeros(num_saidas, num_estados) eye(num_saidas,num_saidas)];
sys_integral = ss(Aaug,Baug,Caug,0,-1);

%Respostas do sistema.
figure(1)
step(sys)
figure(2)
impulse(sys_integral) %Impulso é a derivada do degrau.
