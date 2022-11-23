function [q_next] = propagate_quaternion(q, w, dt)
% Discrete time quaternion propagation 

omega_tol = 1e-5;
n = norm(w);
if n > omega_tol
    c = cos(0.5*n*dt);
    s = sin(0.5*n*dt)/n;
    x = w(1)*s;
    y = w(2)*s;
    z = w(3)*s;
    Omega = [c, z, -y, x; -z, c, x, y; y, -x, c, z; -x, -y, -z, c];
    q_next = Omega*q;
else
    q_next = q;
end
end