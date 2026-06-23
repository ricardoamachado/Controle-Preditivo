function u_limitado = limitar_u(u)
if u > 1
    u_limitado = 1;
elseif u < 0
    u_limitado = 0;
else
    u_limitado = u;
end