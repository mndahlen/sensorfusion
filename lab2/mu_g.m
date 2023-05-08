function [x, P] = mu_g(x, P, yacc, Ra, g0)
    Q = Qq(x);
    dQ = dQqdq(x);
    H = (Q')*g0;
    dH = dQ*g0;

    S = Ra + dH*P*(dH');
    K = P*(dH')*S^(-1);
    eps = yacc - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end