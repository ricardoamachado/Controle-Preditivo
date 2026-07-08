clear;
clc;

%Parâmetros de projeto.
Vs = 48;
Vo = 12;
freq = 50e3;
ripple_corrente = 0.3;
ripple_tensao = 0.02;
Pout = 6;

%Número de espiras.
Npri = 4;
Nsec = 1;
n = Nsec / Npri;

%Duty Cycle, R e corrente no indutor.
D = Vo/(Vs*n + Vo);
R = Vo^2/Pout;
Ilm = (n * Vo) / (R * (1-D));

%Determinação de Lm e C.
C = D / (R * freq * ripple_tensao);
Lm = ((1-D) ^ 2 * R) / (ripple_corrente * n^2 * freq);

%Modelo em espaço de estados.
A_ss = [0 -(1-D)/(n*Lm) ; (1-D)/(n*C) -1/(R*C)];
B_ss = [(Vo/n + Vs)/Lm; -Ilm/(n*C)];
C_ss = [0 1];
D_ss = 0;
sys = ss(A_ss,B_ss,C_ss,D_ss);

%Resposta ao degrau em malha aberta.
step(sys,'b')

% Discretização.
T_sample = 2/freq;
sys_disc = c2d(sys,T_sample,'zoh');
%save("sys_disc.mat","sys_disc")