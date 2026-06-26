function dx = modelo_medio(t, x, u)
% u = Duty cycle
% x(1) = Corrente no indutor.
% x(2) = Tensão no capacitor.
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
    
    %Duty Cycle, R.
    D = Vo/(Vs*n + Vo);
    R = Vo^2/Pout;
  
    %Determinação de Lm e C.
    C = D / (R * freq * ripple_tensao);
    Lm = (R*(1-D) ^ 2) / (ripple_corrente * n^2 * freq);
    dx = zeros(2,1);
    %Derivada da Corrente no indutor.
    dx(1) = -(1-u)/(n*Lm) * x(2) + (u * Vs)/Lm;
    %Derivada da Tensão no capacitor.
    dx(2) = (1-u)/(n*C) * x(1) - x(2)/(R*C);
end