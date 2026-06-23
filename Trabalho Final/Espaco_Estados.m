clc
clear

%% Definindo o sistema em tempo discreto.
load("sys_disc.mat");
A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;
T_sample = sys_disc.Ts;

N_estados = size(A,1);
N_entradas = 1;
N_saidas = 1;
% Aumentando o sistema com a ação integral.
Aaug = [A zeros(N_estados,N_saidas);C*A eye(N_saidas)];
Baug = [B;C*B];
Caug = [zeros(N_saidas, N_estados) eye(N_saidas,N_saidas)];
sys_aug = ss(Aaug,Baug,Caug,0,T_sample);

%% Matrizes de ponderação.
Q = 1*eye(N_saidas);
R = 50*eye(N_entradas);

%% Calculando matrizes de predição.
%Horizontes.
N_predicao = 50;
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

%% Definições das restrições.
% Restrição no sinal de controle (razão cíclica).
u_min = 0;
u_max = 1;
% Inequação na forma F_u * u <= G_u.
F_u = [eye(N_controle);-eye(N_controle)];
G_u = [u_max*ones(N_controle,1);u_min*ones(N_controle,1)];
% Escrevendo a inequação como F_du * delta_u <= G_du.
M_u_para_delta_u = tril(ones(N_controle));
F_du = F_u * M_u_para_delta_u;
% Matriz G_du é calculada dentro da simulação do sistema.

%% Simulação do Sistema Linear.
% Matriz H do quadprog.
H = 2*(B_pred' * C_pred' * Q_pred * C_pred * B_pred + R_pred);
H = (H + H')/2;

%Parâmetros da Simulação
N_sim = 150;
ref = 15;
options = optimoptions('quadprog','Display','off');
W = repmat(ref, N_predicao, 1);
USAR_RESTRICOES = true;

% Armazenar valores da simulação.
y_sim = zeros(1, N_sim);
u_sim = zeros(1, N_sim);
x_sim = zeros(N_estados, N_sim);
delta_u_sim = zeros(1, N_sim);

% Condições Iniciais.
x_sim(:,1) = zeros(N_estados, 1);
u_anterior = 0;
x_anterior = x_sim(:,1);

% Loop da simulação linear.
for k = 1:N_sim
    % Adiciona perturbação na saída a partir da metade da simulação.
    if k > N_sim/2
        perturbacao = 0.5*ref;
    else
        perturbacao = 0;
    end

    % Determinação da saída e do vetor de estados aumentado.
    y_sim(:,k) = C * x_sim(:,k) + perturbacao;
    delta_x = x_sim(:,k) - x_anterior;
    x_aug = [delta_x ; y_sim(:,k)];

    % Matriz F do quadprog.
    F = 2*(C_pred*A_pred*x_aug - W)' * Q_pred *C_pred*B_pred;
    
    % Determinação da matriz G_du das restrições.
    L_du = u_anterior * ones(N_controle, 1);
    G_du = G_u - F_u*L_du;
          
    % Determinação do vetor delta_U.
    if USAR_RESTRICOES
        delta_U = quadprog(H, F, F_du, G_du, [], [], [], [], [], options);
    else
        delta_U = quadprog(H, F, [], [], [], [], [], [], [], options);
    end

    % Proteção contra falha do solver.
    if isempty(delta_U)
        delta_U = zeros(N_controle, 1);
    end
    
    % Armazena o incremento de controle utilizado.
    delta_u_sim(:,k) = delta_U(1:N_entradas);
    
    % Determinação da ação de controle atual.
    if USAR_RESTRICOES
        u_sim(:,k) = u_anterior + delta_u_sim(:,k);
    else
       % Limite físico - Duty cycle entre 0 e 1.
       u_sim(:,k) = limitar_u(u_anterior + delta_u_sim(:,k)); 
    end

    % Simulação da planta.
    if k < N_sim
        x_sim(:,k+1) = A * x_sim(:,k) + B * u_sim(:,k);
    end
    
    % Atualização das memórias para a próxima iteração
    u_anterior = u_sim(:,k);
    x_anterior = x_sim(:,k);
end

%% Simulação Não Linear.

% Ponto de equilibrio do conversor.
U_eq = 0.5;   % Duty cycle de equilíbrio
Y_eq = 12;    % Tensão de saída
X_eq = [0.25; 12]; % Estados no equilíbrio

% 0 <= u_lin + U_eq <= 1
% Limites da razão cíclica em relação ao equilibrio.
u_min_relativo = u_min - U_eq;
u_max_relativo = u_max - U_eq;
G_u = [u_max_relativo*ones(N_controle,1); -u_min_relativo*ones(N_controle,1)];

% Condições Iniciais
x_sim_nl(:,1) = X_eq;
u_anterior_relativo = 0;
x_anterior_relativo = zeros(N_estados, 1);

% Definição da referência em 15 V. Perto do equilibrio.
ref_abs = 15; 
ref_lin = ref_abs - Y_eq;
W = repmat(ref_lin, N_predicao, 1);

% Armazenar valores
y_sim_nl = zeros(1, N_sim);
u_sim_nl = zeros(1, N_sim);
delta_u_sim_nl = zeros(1, N_sim);

% Loop da simulação
for k = 1:N_sim
    
    % Medição da planta e Conversão para Desvio
    y_sim_nl(:,k) = x_sim_nl(2,k);
    % Valor da saída em relação ao equilibrio.
    y_sim_lin = y_sim_nl(:,k) - Y_eq;
    % Determinação do vetor de estados em relação ao equilibrio.
    x_lin_atual = x_sim_nl(:,k) - X_eq; 
    delta_x = x_lin_atual - x_anterior_relativo;
    x_aug = [delta_x ; y_sim_lin];
    
    % Matriz F do quadprog.
    F = 2*(C_pred*A_pred*x_aug - W)' * Q_pred * C_pred * B_pred;
    
    % Determinação da matriz G_du das restrições.
    L_du = u_anterior_relativo * ones(N_controle, 1);
    G_du = G_u - F_u*L_du;
          
    % Determinação do vetor delta_U.
    delta_U = quadprog(H, F, F_du, G_du, [], [], [], [], [], options);
    
    % Caso haja erro no quadprog, o vetor delta_U é nulo.
    if isempty(delta_U)
        delta_U = zeros(N_controle, 1);
    end
    
    delta_u_sim_nl(:,k) = delta_U(1:N_entradas);
    
    % Ação de controle em relação ao equilibrio.
    u_sim_lin_k = u_anterior_relativo + delta_u_sim_nl(:,k);
    u_sim_nl(:,k) = U_eq + u_sim_lin_k; 
    
    % Simulação da planta não linear.
    if k < N_sim
        t_span = [(k-1)*T_sample, k*T_sample];
        [t_ode, x_ode] = ode45(@(t, x) modelo_medio(t, x, u_sim_nl(:,k)), t_span, x_sim_nl(:,k));
        x_sim_nl(:,k+1) = x_ode(end, :)';
    end
    
    % Atualização das memórias para a próxima iteração
    u_anterior_relativo = u_sim_lin_k;
    x_anterior_relativo = x_lin_atual;
end

%% Plotagem dos Resultados
t = (0:N_sim-1) * T_sample;
% Simulação Linear.
figure(1);
subplot(3,1,1);
stairs(t, y_sim, 'k', 'LineWidth', 1.5); hold on;
stairs(t, repmat(ref, 1, N_sim), 'r--', 'LineWidth', 1.5);
title('Tensão de saída vs Referência');
ylabel('Tensão (V)'); grid on; legend('Saída', 'Referência', 'Location', 'best');

subplot(3,1,2);
stairs(t, u_sim, 'k', 'LineWidth', 1.5); hold on;
title('Ação de Controle (u)');
xlabel('Tempo (s)'); ylabel('Duty Cycle'); grid on;
subplot(3,1,3);

stairs(t, delta_u_sim, 'k', 'LineWidth', 1.5); hold on;
title('Incremento de controle (Δu)');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on;

% Simulação não Linear.
figure(2);
subplot(3,1,1);
stairs(t, y_sim_nl, 'k', 'LineWidth', 1.5); hold on;
stairs(t, repmat(ref_abs, 1, N_sim), 'r--', 'LineWidth', 1.5);
title('Tensão de saída vs Referência');
ylabel('Tensão (V)'); grid on; legend('Saída', 'Referência', 'Location', 'best');

subplot(3,1,2);
stairs(t, u_sim_nl, 'k', 'LineWidth', 1.5); hold on;
title('Ação de Controle (u)');
xlabel('Tempo (s)'); ylabel('Duty Cycle'); grid on;
subplot(3,1,3);

stairs(t, delta_u_sim_nl, 'k', 'LineWidth', 1.5); hold on;
title('Incremento de controle (Δu)');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on;