clc;
clear;
load("sys_disc.mat")

%% Parâmetros comuns.
N_sim        = 400;
N_truncamento = 150;
T_sample     = sys_disc.Ts;
t            = (0:N_sim-1) * T_sample;
Q            = 1;
R            = 1000;
USAR_RESTRICOES = true;
ref          = 12;

controladores    = {'GPC', 'DMC', 'SS-MPC'};
n_ctrl           = numel(controladores);
ylabels_subplots = {'Tensão (V)', 'Duty Cycle', 'Incremento de controle (Δu)'};


%% Variação de R.
N_predicao_R = 50;
N_controle_R = 50;
R_valores    = [10, 100, 500, 1000, 5000];
legend_R     = arrayfun(@(r) sprintf('R = %d', r), R_valores, 'UniformOutput', false);

run_fns_R = { ...
    @run_gpc_linear, ...
    @run_dmc_linear, ...
    @run_ssmpc_linear };

build_args_R = @(R_val) { ...
    {sys_disc, ref, N_predicao_R, N_controle_R, N_sim, Q, R_val,  USAR_RESTRICOES}, ...
    {sys_disc, ref, N_predicao_R, N_controle_R, N_sim, N_truncamento, Q, R_val, USAR_RESTRICOES}, ...
    {sys_disc, ref, N_predicao_R, N_controle_R, N_sim, Q, R_val,  USAR_RESTRICOES} };

plot_variacao(run_fns_R, R_valores, legend_R, 1, build_args_R, ...
              controladores, ylabels_subplots, n_ctrl, N_sim, t, ref);

%% Variação de N.
N_valores   = [2, 5, 10, 20, 40];
legend_N    = arrayfun(@(n) sprintf('N = %d', n), N_valores, 'UniformOutput', false);

run_fns_N = { ...
    @run_gpc_linear, ...
    @run_dmc_linear, ...
    @run_ssmpc_linear };

build_args_N = @(N_val) { ...
    {sys_disc, ref, N_val, N_val, N_sim, Q, R, USAR_RESTRICOES}, ...
    {sys_disc, ref, N_val, N_val, N_sim, N_truncamento, Q, R, USAR_RESTRICOES}, ...
    {sys_disc, ref, N_val, N_val, N_sim, Q, R, USAR_RESTRICOES} };

plot_variacao(run_fns_N, N_valores, legend_N, 4, build_args_N, ...
              controladores, ylabels_subplots, n_ctrl, N_sim, t, ref);

%% Variação do horizonte de truncamento do DMC.
N_predicao_T  = 50;
N_controle_T  = 50;
NT_valores    = [20, 50, 100, 150];
legend_NT     = arrayfun(@(n) sprintf('N_{trunc} = %d', n), NT_valores, 'UniformOutput', false);

n_NT = numel(NT_valores);

% Pré-alocação.
y_nt       = zeros(N_sim, n_NT);
u_nt       = zeros(N_sim, n_NT);
delta_u_nt = zeros(N_sim, n_NT);

% Simulação — apenas DMC varia.
for k = 1:n_NT
    [y_nt(:,k), u_nt(:,k), delta_u_nt(:,k)] = ...
        run_dmc_linear(sys_disc, ref, N_predicao_T, N_controle_T, ...
                       N_sim, NT_valores(k), Q, R, USAR_RESTRICOES);
end

% Plot — figura única para o DMC.
figure(7);
clf;
sgtitle('DMC — Variação de N_{truncamento}', 'FontWeight', 'bold', 'FontSize', 14);
dados_nt = {y_nt, u_nt, delta_u_nt};

for s = 1:3
    subplot(3, 1, s);
    hold on;
    for k = 1:n_NT
        stairs(t, dados_nt{s}(:,k), 'LineWidth', 1.8);
    end
    if s == 1
        yline(ref, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Referência');
    end
    hold off;
    ylabel(ylabels_subplots{s});
    if s == 3
        xlabel('Tempo (s)');
    end
    legend(legend_NT, 'Location', 'best');
    grid on;
    box on;
end

%% Variação de N_predicao.
N_controle = [2, 5, 10, 20, 40];
N_predicao = 50;
legend_N     = arrayfun(@(n) sprintf('N_u = %d', n), N_controle, 'UniformOutput', false);

run_fns_N = { ...
    @run_gpc_linear, ...
    @run_dmc_linear, ...
    @run_ssmpc_linear };

build_args_N = @(N_ctrl_val) { ...
    {sys_disc, ref, N_predicao, N_ctrl_val, N_sim, Q, R,  USAR_RESTRICOES}, ...
    {sys_disc, ref, N_predicao, N_ctrl_val, N_sim, N_truncamento, Q, R, USAR_RESTRICOES}, ...
    {sys_disc, ref, N_predicao, N_ctrl_val, N_sim, Q, R,  USAR_RESTRICOES} };

plot_variacao(run_fns_N, N_controle, legend_N, 8, build_args_N, ...
              controladores, ylabels_subplots, n_ctrl, N_sim, t, ref);

%% Função auxiliar para plotagem dos resultados.

    function plot_variacao(run_fns, param_valores, legend_str, ...
                           fig_offset, build_args_fn, ...
                           controladores, ylabels_subplots, ...
                           n_ctrl, N_sim, t, ref)

        n_param = numel(param_valores);

        % Pré-alocação.
        y       = cell(n_ctrl, 1);
        u       = cell(n_ctrl, 1);
        delta_u = cell(n_ctrl, 1);
        for c = 1:n_ctrl
            y{c}       = zeros(N_sim, n_param);
            u{c}       = zeros(N_sim, n_param);
            delta_u{c} = zeros(N_sim, n_param);
        end

        % Simulação.
        for k = 1:n_param
            args = build_args_fn(param_valores(k));   % cell array de args
            for c = 1:n_ctrl
                [yc, uc, duc] = run_fns{c}(args{c}{:});
                y{c}(:,k)       = yc;
                u{c}(:,k)       = uc;
                delta_u{c}(:,k) = duc;
            end
        end

        % Plot.
        for c = 1:n_ctrl
            figure(fig_offset + c - 1);
            clf;
            sgtitle(controladores{c}, 'FontWeight', 'bold', 'FontSize', 14);
            dados = {y{c}, u{c}, delta_u{c}};

            for s = 1:3
                subplot(3, 1, s);
                hold on;
                for k = 1:n_param
                    stairs(t, dados{s}(:,k), 'LineWidth', 1.8);
                end
                if s == 1
                    yline(ref, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Referência');
                end
                hold off;
                ylabel(ylabels_subplots{s});
                if s == 3
                    xlabel('Tempo (s)');
                end
                legend(legend_str, 'Location', 'best');
                grid on;
                box on;
            end
        end
    end
