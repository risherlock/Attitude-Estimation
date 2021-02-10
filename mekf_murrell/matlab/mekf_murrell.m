function [q_hat, b_hat, dx_hat, P_k] = ...
    mekf_murrell(q, b, w, P, mb, mr, sigma_r, sigma_g, sigma_b, dt)
%%% Multiplicative Extended Kalman Filter with Murrell's algorithm.
%
% Inputs:
%   q   = Quaternion estimate of previous time step [q1,q2,q3,q4]' [4x1]
%   b   = Gyro-bias estimate of previous time step [3x1]
%   w   = Updated angular velocity (i.e. w = gyro_w - b) [3x1]
%   P   = State covariance update of previous time step [6x6]
%   mb  = n sensor measurement vectors in body frame [3xn]
%   mr  = n measurement vectors in inertial frame [3xn]
%   sigma_r = Vector sensors' noise standard deviations [1xn]
%   sigma_g = Gyro noise standard deviation matrix, rad/sec^(1/2) [3x3]
%   sigma_b = Gyro-bias noise standard deviation matrix, rad/sec^(3/2) [3x3]
%   Q   = Process noise covariance [6x6]
%   dt  = Sample time (sec)
%
% Outputs:
%   q_hat   = Updated quaternion estimate [4x1]
%   b_hat   = Updated gyro-biases estimate [3x1]
%   dx_hat  = Updated state estimate [6x1]
%   P_k     = Updated estimate covariance [6x6]
%
% References:
%   [1] Markley, Crassidis - Fundamentals of Spacecraft Attitude
%       Determination and Control (2014)
%   [2] Crassidis, Junkins - Optimal Estimation of Dynamic Systems
%        (2nd ed.) (2011)
%
% Note:
%    1. The equations mentioned in the comments references to [2].
%    2. For 3-axis gyro, sigma_g = diag[sigma_gx, sigma_gy, sigma_gz] and 
%       sigma_b = diag[sigma_bx, sigma_by, sigma_bz].
%
% Rishav (2021/1/18)

% Propagate state, gyro-bias, quaternion and estimate covariance
% Note: Gyro bias propagation is constant. Eq.(7.42b)(i.e. b = b)
Q = computeProcessCovariance(sigma_g, sigma_b, dt); % Eq.(7.46)
q = propQuaternion(q, w, dt); % Eq.(7.39)
P = propStateCovariance(P, Q, w, dt); % Eq.(7.43)
dx = [0,0,0,0,0,0]';
A = quaternion2A(q); % Attitude matrix

% Murrell's loop
[~,n] = size(mb);
for i_iters = 1 : n
    % Sensitivity matrix (H[3x6]) and residal (e[3x1])
    Ar  = A*mr(:,i_iters); % m_b = A * m_r (ideally)
    ArX = [0, -Ar(3), Ar(2); Ar(3), 0, -Ar(1); -Ar(2),  Ar(1), 0];
    H   = [ArX, zeros(3)]; % [3x6] Eq.(7.29)
    e   = mb(:,i_iters) - Ar;
    
    % MEKF update
    K = P*H'/(H*P*H' + sigma_r(i_iters)^2 * eye(3)); % [6x3]
    dx = dx + K*(e - H*dx); % Eq.(7.30)
    
    P = (eye(6) - K*H)*P;
end
dx_hat = dx;
P_k = P;

% Update q and b
dq = dx_hat(1:3);
db = dx_hat(4:6);
q_hat = upQuaternion(q, dq); % Eq.(7.34)
b_hat = b + db; % Eq.(7.32)
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Quaternion to rotation matrix
function [A] = quaternion2A(q)
A = zeros(3);
A(1,1) = + q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
A(2,2) = - q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2;
A(3,3) = - q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;
A(1,2) = 2*(q(1)*q(2) + q(3)*q(4));
A(1,3) = 2*(q(1)*q(3) - q(2)*q(4));
A(2,1) = 2*(q(1)*q(2) - q(3)*q(4));
A(2,3) = 2*(q(2)*q(3) + q(1)*q(4));
A(3,1) = 2*(q(1)*q(3) + q(2)*q(4));
A(3,2) = 2*(q(2)*q(3) - q(1)*q(4));
end

% Discrete-time quaternion propagation
function [q_out] = propQuaternion(q, w, dt)
omega_tol = 1e-5;
n = norm(w);
if n > omega_tol
    % Eq.(7.40)
    c = cos(0.5*n*dt);
    s = sin(0.5*n*dt)/n;
    x = w(1)*s;
    y = w(2)*s;
    z = w(3)*s;
    Omega = [c, z, -y, x; -z, c, x, y; y, -x, c, z; -x, -y, -z, c];
    q_out = Omega*q; % Eq.(7.39)
else
    q_out = q;
end
end

% Update quaternion using error-quaternion
function [q_out] = upQuaternion(q, del_alpha)
Xi = [q(4), -q(3), q(2); ...
    q(3), q(4), -q(1); ...
    -q(2), q(1), q(4); ...
    -q(1), -q(2), -q(3)];
q_out = q + 0.5*Xi*del_alpha; % Eq.(7.34)
q_out = q_out/norm(q_out); % Quaternion normalization
end

% Discrete propagation of the state error covariance
function [P_out] = propStateCovariance(P, Q, w, dt)
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

% Q using gyro and gyro-bias noise standard deviations
function [Q] = computeProcessCovariance(sigma_g, sigma_b, dt)
% Inputs:
%   sigma_g = Gyro noise standard deviation, rad/sec^(1/2) [3x3]
%   sigma_b = Gyro-bias noise standard deviation, rad/sec^(3/2) [3x3]

% Eq.(7.46)
Q = zeros(6); % Discrete Q
ssg = sigma_g.^2;
ssb = sigma_b.^2;
Q(1:3,1:3) = ssg*dt + ssb*dt*dt*dt/3;
Q(1:3,4:6) = 0.5*ssb*dt*dt;
Q(4:6,1:3) = Q(1:3,4:6)';
Q(4:6,4:6) = dt*ssb;
end
