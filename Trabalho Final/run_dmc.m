function [y_sim, u_sim, delta_u_sim, y_sim_nl, u_sim_nl, delta_u_sim_nl] = run_dmc(sys_disc,ref_abs,ref,N_predicao,N_controle,N_sim,N_truncamento,Q,R,USAR_RESTRICOES)


A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;
D = sys_disc.D;
T_sample = sys_disc.Ts;

% Horizontes.
N_estados = size(A,1);
N_amostras = N_predicao + N_truncamento;

[g_vetor, t_vetor] = step(sys_disc, 0:T_sample:N_amostras*T_sample);
% Corta o elemento associado com k = 0.
g_vetor(1) = [];
t_vetor(1) = [];

G = zeros(N_predicao, N_controle);
G(:, 1) = g_vetor(1:N_predicao);

for i = 1:N_controle-1
    % Desloca g_vetor i unidades e preenche início com zeros.
    g_vetor_i = circshift(g_vetor, i);
    g_vetor_i(1:i) = 0;
    G(:, i+1) = g_vetor_i(1:N_predicao);  % Atribui coluna já recortada
end

% Monta a matriz da resposta livre Gf.
Gf = zeros(N_predicao, N_truncamento);

for passo = 1:N_predicao
    gf_linha = diff_n(g_vetor, passo)';
    Gf(passo, :) = gf_linha(1:N_truncamento);
end

%% Matrizes de ponderação.
Q_pred = Q * eye(N_predicao);
R_pred = R * eye(N_controle);
H_qp = 2*(R_pred + G'*Q_pred*G);
H_qp = (H_qp + H_qp')/2;

%% Definições das restrições.
% Restrição no sinal de controle (razão cíclica).
u_min = 0;
u_max = 1;
% Inequação na forma F_u * u <= G_u.
F_u = [eye(N_controle);-eye(N_controle)];
G_u = [u_max*ones(N_controle,1);-u_min*ones(N_controle,1)];
% Escrevendo a inequação como F_du * delta_u <= G_du.
M_u_para_delta_u = tril(ones(N_controle));
F_du = F_u * M_u_para_delta_u;
% Matriz G_du é calculada dentro da simulação do sistema.

%% Simulação Linear.
% Vetores utilizados no cálculo da resposta livre.
delta_u_passado = zeros(N_truncamento,1);

%Parâmetros da Simulação
options = optimoptions('quadprog','Display','off');
W = repmat(ref, N_predicao, 1);

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
        perturbacao = 0.2*ref;
    else
        perturbacao = 0;
    end

    % Determinação da saída.
    y_sim(:,k) = C * x_sim(:,k) + perturbacao;

    % Vetor da resposta livre.
    f_livre = y_sim(:,k) + Gf * delta_u_passado;

    % Matriz F do quadprog.
    F_qp = 2*(f_livre - W)' * Q_pred * G;
    F_qp = F_qp';
    % Determinação da matriz G_du das restrições.
    L_du = u_anterior * ones(N_controle, 1);
    G_du = G_u - F_u*L_du;
          
    % Determinação do vetor delta_U.
    if USAR_RESTRICOES
        delta_U = quadprog(H_qp, F_qp, F_du, G_du, [], [], [], [], [], options);
    else
        delta_U = quadprog(H_qp, F_qp, [], [], [], [], [], [], [], options);
    end

    % Proteção contra falha do solver.
    if isempty(delta_U)
        delta_U = zeros(N_controle, 1);
    end
     
    % Armazena o incremento de controle utilizado.
    delta_u_sim(:,k) = delta_U(1);
    delta_u_passado = circshift(delta_u_passado,1);
    delta_u_passado(1) = delta_u_sim(:,k);

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
delta_u_passado_nl = zeros(N_truncamento,1);


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

% Definição da referência não linear. Perto do equilibrio.
ref_lin = ref_abs - Y_eq;
W = repmat(ref_lin, N_predicao, 1);

% Armazenar valores
y_sim_nl = zeros(1, N_sim);
u_sim_nl = zeros(1, N_sim);
delta_u_sim_nl = zeros(1, N_sim);

% Loop da simulação
for k = 1:N_sim
    if k > N_sim * 0.6
        perturbacao = 1;
    else
        perturbacao = 0;
    end
    % Medição da planta e Conversão para Desvio
    y_sim_nl(:,k) = x_sim_nl(2,k) + perturbacao;

    % Valor da saída em relação ao equilibrio.
    y_sim_relativo = y_sim_nl(:,k) - Y_eq;

    % Determinação do vetor de estados em relação ao equilibrio.
    x_lin_atual = x_sim_nl(:,k) - X_eq; 

    % Vetor da resposta livre.
    f_livre = y_sim_relativo + Gf * delta_u_passado_nl;

    % Matriz F do quadprog.
    F_qp = 2*(f_livre - W)' * Q_pred * G;
    F_qp = F_qp';
    
    % Determinação da matriz G_du das restrições.
    L_du = u_anterior_relativo * ones(N_controle, 1);
    G_du = G_u - F_u*L_du;
          
    % Determinação do vetor delta_U.
    delta_U = quadprog(H_qp, F_qp, F_du, G_du, [], [], [], [], [], options);
    
    % Caso haja erro no quadprog, o vetor delta_U é nulo.
    if isempty(delta_U)
        delta_U = zeros(N_controle, 1);
    end
    
    % Armazena o incremento de controle utilizado.
    delta_u_sim_nl(:,k) = delta_U(1);
    delta_u_passado_nl = circshift(delta_u_passado_nl,1);
    delta_u_passado_nl(1) = delta_u_sim_nl(:,k);
    
    % Ação de controle em relação ao equilibrio.
    u_sim_relativo_k = u_anterior_relativo + delta_u_sim_nl(:,k);
    u_sim_nl(:,k) = U_eq + u_sim_relativo_k; 
    
    % Simulação da planta não linear.
    if k < N_sim
        t_span = [(k-1)*T_sample, k*T_sample];
        [t_ode, x_ode] = ode45(@(t, x) modelo_medio(t, x, u_sim_nl(:,k)), t_span, x_sim_nl(:,k));
        x_sim_nl(:,k+1) = x_ode(end, :)';
    end
    
    % Atualização das memórias para a próxima iteração
    u_anterior_relativo = u_sim_relativo_k;
end
end
