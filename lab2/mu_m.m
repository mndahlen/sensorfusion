function [x, P] = mu_m(x, P, ymag, Rm, m0)
    Q = Qq(x);
    [dQ0, dQ1, dQ2, dQ3] = dQqdq(x);
    H = (Q')*m0; 
    dH = [dQ0'*m0, dQ1'*m0, dQ2'*m0, dQ3'*m0];

    S = Rm + dH*P*(dH');
    K = P*(dH')*S^(-1); % (S\eye(length(S))
    eps = ymag - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end