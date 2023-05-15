function [x, P] = mu_g(x, P, yacc, Ra, g0)
    Q = Qq(x);
    [dQ0, dQ1, dQ2, dQ3] = dQqdq(x);
    H = (Q')*g0; 
    dH = [dQ0'*g0, dQ1'*g0, dQ2'*g0, dQ3'*g0]; % Seems to drift?

    S = Ra + dH*P*(dH');
    K = P*(dH')*S^(-1);
    eps = yacc - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end