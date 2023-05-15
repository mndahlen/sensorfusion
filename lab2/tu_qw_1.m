function [x, P] = tu_qw(x, P, omega, T, Rw)
    % If no omega, use previous omega?
    S_q = Sq(x);
    S_omega = Somega(omega);
    I = eye(size(S_omega));

    F = (I + 0.5*S_omega*T);
    G = 0.5*T*S_q;
    Q = Rw;
    
    x = F*x;
    P = F*P*(F') + G*Q*(G');
end

