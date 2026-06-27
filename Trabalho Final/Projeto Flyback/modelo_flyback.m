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
C = D/(R*freq*ripple_tensao);
Lm = R*(1-D)^2/(ripple_corrente * n^2 * freq);

%Modelo em espaço de estados.
A_ss = [0 -(1-D)/(n*Lm) ; (1-D)/(n*C) -1/(R*C)];
B_ss = [Vo/(n*Lm) + Vs/Lm; -Ilm/(n*C)];
C_ss = [0 1];
sys = ss(A_ss,B_ss,C_ss,0);
T_sample = 2/freq; % Ajustar posteriormente.
sys_disc = c2d(sys,T_sample,'zoh');

%Comparação entre simulação e modelo.
load("Resp Freq 0.01 Amplitude.mat")
[mag, phase, freqs] = bode(sys, 2*pi*Frequency);
mag = squeeze(mag);
mag = 20*log10(mag);
phase = squeeze(phase);
phase = phase - 360;
figure(1)
tiledlayout(2,1)
nexttile
semilogx(Frequency,amp_Vo2,'b--','LineWidth', 2.5)
hold on
semilogx(Frequency,mag,'k','LineWidth', 1)
legend('Simulação','Modelo')
xlabel("Frequência (Hz)")
ylabel("Magnitude (dB)")
nexttile
semilogx(Frequency,phase_Vo2,'b--','LineWidth', 2.5)
hold on
semilogx(Frequency,phase,'k','LineWidth', 1)
legend('Simulação','Modelo')
xlabel("Frequência (Hz)")
ylabel("Fase (graus)")