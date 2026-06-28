clc;
clear;

load("sys_disc.mat")

% Horizontes.
N_predicao = 50;
N_controle = 50;
N_sim = 200;
N_truncamento = 150;
% Intervalo de tempo da simulação.
T_sample = sys_disc.Ts;
t = (0:N_sim-1) * T_sample;
% Matrizes de ponderação.
Q = 1;
R = 10000;
USAR_RESTRICOES = true;
% Referência para o modelo não linear.
ref_abs = 14;
% Referência para o modelo linear.
ref = 12;

% Execução dos contraladores.
[y_sim_gpc, u_sim_gpc, delta_u_sim_gpc, y_sim_nl_gpc, u_sim_nl_gpc, delta_u_sim_nl_gpc] = run_gpc(sys_disc,ref_abs,ref,N_predicao,N_controle,N_sim,Q,R,USAR_RESTRICOES);
[y_sim_dmc, u_sim_dmc, delta_u_sim_dmc, y_sim_nl_dmc, u_sim_nl_dmc, delta_u_sim_nl_dmc] = run_dmc(sys_disc,ref_abs,ref,N_predicao,N_controle,N_sim,N_truncamento,Q,R,USAR_RESTRICOES);
[y_sim_mpc, u_sim_mpc, delta_u_sim_mpc, y_sim_nl_mpc, u_sim_nl_mpc, delta_u_sim_nl_mpc] = run_ssmpc(sys_disc,ref_abs,ref,N_predicao,N_controle,N_sim,Q,R,USAR_RESTRICOES);

% Simulação Linear.
figure(1);
subplot(3,1,1);
stairs(t, y_sim_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, y_sim_dmc, 'b', 'LineWidth', 2);
stairs(t, y_sim_mpc, 'r--', 'LineWidth', 3);
title('Tensão de saída');
ylabel('Tensão (V)');
grid on; 
legend('GPC','DMC', 'SSMPC','Location', 'best');

subplot(3,1,2);
stairs(t, u_sim_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, u_sim_dmc, 'b', 'LineWidth', 2);
stairs(t, u_sim_mpc, 'r--', 'LineWidth', 2);
legend('GPC','DMC', 'SSMPC','Location', 'best')
title('Ação de Controle (u)');
xlabel('Tempo (s)');
ylabel('Duty Cycle');
grid on;

subplot(3,1,3);
stairs(t, delta_u_sim_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, delta_u_sim_dmc, 'b', 'LineWidth', 2);
stairs(t, delta_u_sim_mpc, 'r--', 'LineWidth', 3);
legend('GPC','DMC', 'SSMPC','Location', 'best')
title('Incremento de controle (Δu)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;

% Simulação não Linear.
figure(2);
subplot(3,1,1);
stairs(t, y_sim_nl_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, y_sim_nl_dmc, 'b', 'LineWidth', 2);
stairs(t, y_sim_nl_mpc, 'r--', 'LineWidth', 3);
title('Tensão de saída');
ylabel('Tensão (V)');
grid on; 
legend('GPC','DMC', 'SSMPC','Location', 'best');

subplot(3,1,2);
stairs(t, u_sim_nl_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, u_sim_nl_dmc, 'b', 'LineWidth', 2);
stairs(t, u_sim_nl_mpc, 'r--', 'LineWidth', 3);
legend('GPC','DMC', 'SSMPC','Location', 'best')
title('Ação de Controle (u)');
xlabel('Tempo (s)');
ylabel('Duty Cycle');
grid on;

subplot(3,1,3);
stairs(t, delta_u_sim_nl_gpc, 'k', 'LineWidth', 2);
hold on;
stairs(t, delta_u_sim_nl_dmc, 'b', 'LineWidth', 2);
stairs(t, delta_u_sim_nl_mpc, 'r--', 'LineWidth', 3)
legend('GPC','DMC', 'SSMPC','Location', 'best')
title('Incremento de controle (Δu)');
xlabel('Tempo (s)');
ylabel('Amplitude');
grid on;