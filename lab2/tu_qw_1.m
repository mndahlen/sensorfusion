function [x, P] = tu_qw(x, P, omega, T, Rw)
% If no omega, use previous omega?
    I = eye(size(Rw));
    S_q = Sq(x);
    S_omega = Somega(omega);
    
    F_k = (I + 0.5*S_omega*T);
    G_k = 0.5*T*S_q;
    Q = Rw;

    x = F*x;
    P = F*P*(F_k') + G_k*Q*(G_k');
end

