function [x_hat,P_hat] = UKF(x_prev,P_prev,z,dt,sys_params,ukf_params)
% Extended Kalman filter implementation for noisy pendulum tracking
%
% Input: n =
% x = Present state i.e. previous estimate (nX1)
% u = Control input (px1)
% z = Measurement (mx1)
% sys_params = {F,B,H,P,Q,R}, a cell array
%       F = State transition matrix (nXn)
%       B = Control matrix that maps control to state variables (nXp)
%       H = Measurement matrix that maps measurements to state (mXn)
%       P = State covariance matrix (n*n)
%       Q = Process covariance matrix (nX1)
%       R = Measurement covariance matrix (mXm)
%
% Output:
%   state       = State estimate by KF (nX1)
%   covariance  = State covariance estimate by KF (nXn)
%
% Rishav (2020/9/3)

% Unpack system parameters
g = sys_params{1};
L = sys_params{2};
R = sys_params{3};
Q = sys_params{4};
d = sys_params{5};

n = 2; % Size of the state vector
alpha  = 1; % Primary scaling parameter
beta = 2; % Secondary scaline parameter
kappa = 0; % Tertiary scaling parameter

lambda = alpha^2*(n+kappa) - n;
W_mean = ones(2*n+1,1)*1/(2*(n+lambda));
W_cov = W_mean; W_mean(1) = lambda/(lambda+n);
W_cov(1) = lambda/(lambda+n) + 1 - alpha^2 + beta;

m = 3;

% Prediction
sqrtP = chol(P_prev,'lower');
chi = [x_prev, ... % 2n+1 sigma points
    x_prev + sqrt(n+lambda)*sqrtP, ...
    x_prev - sqrt(n+lambda)*sqrtP];

% f = 2(a + 1) (Eqn 21)
a = 1;
f = 2*(a + 1);

chi_x = zeros(n,2*n+1); % Propagated sigma points
chi_z = zeros(m,2*n+1);
chi_x(1) = [delp_prev', beta_prev']'; % Eqn 5b

% Propagate sigma points (chi(0))
% for chi(:, 1)
q = q_prev;
q0_

% Propagate sigma points (chi(1) to chi(12))
for i_iters = 2:2*n+1
    % Extract states from sigma points
    del_p = chi(1:3, i_iters);
    beta  = chi(4:6, i_iters);  
    
    % Compute del_q from del_p
    del_q4  = (-a * norm(del_p)^2 + f*sqrt(f^2 + (1 - a^2)*del_p^2))/ ...
              f^2 + norm(del_p)^2; % Eqn 21a and 33a
    del_rho = (1/f)*(a + del_q4)*del_p; % Eqn 21 and eqn 33 
    del_q   = [del_rho', del_q4];
    
    % Compute new q using q_prev and del_q
    q = quat_prod(del_q, q_prev); % Eqn 32b
    
    % Propogate quaternion
    omega   = omega_gyro - beta; % Eqn 35 and Eqn 25a
    norm_w  = norm(omega);
    cosine  = cos(0.5*norm_w*dt);  
    psi     = sin(0.5*norm_w*dt)*omega/norm_w;
    psiX    = [0, -psi(3), psi(2); psi(3), 0, -psi(1); -psi(2), psi(1), 0];
    mat     = [cosine*eye(3) - psiX,  psi; -psi', cosine]; % Eqn 29
    q_      = mat*q; % Eqn 34 and Eqn 28
    
    % Compute propagated error quaternion
    del_q_ = quat_prod(q_, q0_);
    
    % Compute propagated del_p
    del_p_ = f*del_q_(1:3)/(a + del_q_(4));
    
    chi_x(1:3, i_iters) = del_p_;
    chi_x(4:6, i_iters) = beta;
    chi_z(:,i_iters) = sensorModel();
end
x_ = chi_x*W_mean; % Predicted state
z_ = chi_z*W_mean; % Predicted measurement

P_ = Q; Pyy = R; Pxy = zeros(n,m);
for i_iters = 1:2*n+1
    M = chi_x(:,i_iters) - x_;
    N = chi_z(:,i_iters) - z_;
    
    P_ = P_ + W_cov(i_iters)*(M*M'); % Predicted covariance
    Pyy = Pyy + W_cov(i_iters)*(N*N'); % Measurement/innovation covariance
    Pxy = Pxy + W_cov(i_iters)*(M*N'); % State-measurement cross-covariance
end

%%% Update
K = Pxy/Pyy; % Kalman gain
y_pre = z - z_; % Innovation or measurement pre-fit residual
x_hat = x_ + K*y_pre; % Updated state estimate
P_hat = P_ - K*Pyy*K'; % Updated estimate covariance
end

% Propagate GRP error dynamics
function [a] = dynamicsModel()
a = 0;
end

% Propagate gyro model
function [a] = sensorModel()
a = 0;
end



