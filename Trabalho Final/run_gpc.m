function [y_sim, u_sim, delta_u_sim, y_sim_nl, u_sim_nl, delta_u_sim_nl] = run_gpc(sys_disc,ref_abs,ref,N_predicao,N_controle,N_sim,Q,R,USAR_RESTRICOES)

A = sys_disc.A;
B = sys_disc.B;
C = sys_disc.C;
D = sys_disc.D;
T_sample = sys_disc.Ts;
N_estados = size(A,1);

[B_poly, A_poly] = ss2tf(A,B,C,D);
% Função de transf. na forma (z^-1)*(B(z^-1))/(A(z^-1))
B_poly = B_poly(2:end); % Retira o z^-1.

% Multiplica A por (1-z^-1) para obter o A~.
A_til = conv(A_poly, [1, -1]);

%Cada linha j de E representa o polinômio E_j.
%Cada linha j de F representa o polinômio F_j.
[E, F] = Diofantina(A_til,1,N_predicao);

% Calcula o tamanho de cada linha após a convolução.
tamanho_E = size(E, 2);
tamanho_B = length(B_poly);
tamanho_conv = tamanho_E + tamanho_B - 1;

M = zeros(N_predicao, tamanho_conv);
% Matriz H multiplica os incrementos de controle passados.
% N linhas e n_b colunas, com n_b grau do polinômio B(z^-1).
H = zeros(N_predicao,tamanho_B-1);
%Calcula o produto M_j = E_j * B_j e armazena em M.
for idx = 1:N_predicao
    M(idx, :) = conv(E(idx, :), B_poly);
    H(idx,:) = M(idx, idx+1:idx+tamanho_B-1);
end

% Separa o vetor com coeficientes da resp. ao degrau.
g_vetor = M(N_predicao,1:N_predicao);
g_vetor = g_vetor';

% Monta a matriz da resposta forçada G.
G = g_vetor;
for i=1:N_controle-1
    %Desloca o g_vetor i unidades e preenche início com zeros.
    g_vetor_i = circshift(g_vetor,i);
    g_vetor_i(1:i) = 0;
    G = cat(2,G,g_vetor_i);
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
% ATENÇÃO - Primeiro elemento do y_passado é o y(k).
% ATENÇÃO - Primeiro elemento do delta_u_passado é o delta_u(k-1).

y_passado = zeros(size(F,2),1);
delta_u_passado = zeros(size(H,2),1);

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
    y_passado = circshift(y_passado,1);
    y_passado(1) = y_sim(:,k);
    
    % Vetor da resposta livre.
    f_livre = F * y_passado + H * delta_u_passado;

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

% Vetores utilizados no cálculo da resposta livre.
y_passado_nl = zeros(size(F,2),1);
delta_u_passado_nl = zeros(size(H,2),1);

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

% Definição da referência em 14 V. Perto do equilibrio.
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
    y_passado_nl = circshift(y_passado_nl,1);
    y_passado_nl(1) = y_sim_relativo;

    % Determinação do vetor de estados em relação ao equilibrio.
    x_lin_atual = x_sim_nl(:,k) - X_eq; 

    % Vetor da resposta livre.
    f_livre = F* y_passado_nl + H * delta_u_passado_nl;

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
        [~, x_ode] = ode45(@(t, x) modelo_medio(t, x, u_sim_nl(:,k)), t_span, x_sim_nl(:,k));
        x_sim_nl(:,k+1) = x_ode(end, :)';
    end
    
    % Atualização das memórias para a próxima iteração
    u_anterior_relativo = u_sim_relativo_k;
end

end