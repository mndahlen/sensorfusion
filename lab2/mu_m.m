function [x, P] = mu_m(x, P, ymag, Rm, m0)
% correct?
    Q = Qq(x);
    dQ = dQqdq(x);
    H = (Q')*m0;
    dH = dQ*m0;

    S = Rm + dH*P*(dH');
    K = P*(dH')*S^(-1);
    eps = ymag - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end