function [P_out] = propagate_state_covariance(P, Q, w, dt)
% Discrete propagation of the state error covariance

omega_tol = 1e-5;
n = norm(w);

if n > omega_tol
    % Eq.(7.45)
    s = sin(n*dt);
    c = cos(n*dt);
    wX  = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
    Phi_11 = eye(3) - wX*s/n + wX*wX*(1 - c)/n^2;
    Phi_12 = wX*(1 - c)/n^2 - eye(3)*dt - wX*wX*(n*dt - s)/n^3;
else
    % Steady-state Phi
    % Eq.(7.45) and Eq.(7.51b)
    Phi_11 = eye(3);
    Phi_12 = - eye(3)*dt;
end
    Phi_21 = zeros(3);
    Phi_22 = eye(3);
    Phi = [Phi_11, Phi_12; Phi_21, Phi_22];
    
    Y = [-eye(3), zeros(3); zeros(3), eye(3)]; % Eq.(7.23b)
    P_out = Phi*P*Phi' + Y*Q*Y'; % Eq.(7.43)
end
 