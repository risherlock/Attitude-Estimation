% Simulation demonstrating attitude estimation of rigid body using MEKF
% with Murrell's algorithm.
%
% References:
%   [1] Markley, Crassidis - Fundamentals of Spacecraft Attitude
%       Determination and Control (2014)
%   [2] Crassidis, Junkins - Optimal Estimation of Dynamic Systems
%        (2nd ed.) (2011)
%
% Note:
%    1. The equations mentioned in the comments references to [2].
%
% Rishav (2020/1/18)
clc, clear, close all

% Simulation params
sim_time = 1000; % Seconds
dt = 0.01;
time = 0:dt:sim_time;

% Angular velocities
w1_true = 0.1 + 0.1*sin(time+1); % rad/s
w2_true = 0.2 + 0.3*sin(time+2); % rad/s
w3_true = 0.3 + 0.9*sin(time+3); % rad/s

% Initial conditions
q0  = [0,0,0,1]'; % Quaternion
P_a = diag([0.005,0.005,0.005]); % Attitude error covariance, deg^2
P_g = diag([0.005,0.005,0.005]); % Gyro covariance, (deg/hr)^2
dx0 = [0,0,0,0,0,0]'; % State
b0  = [0,0,0]'; % Bias, deg/hr

% Noise standard deviations of 3-axis gyro
sigma_bx = 1e-2; % Gyro-bias noise, rad/sec^(3/2)
sigma_by = sigma_bx;
sigma_bz = sigma_bx;
sigma_gx = 1e-3; % Gyro noise, rad/sec^(1/2)
sigma_gy = sigma_gx;
sigma_gz = sigma_gx;

% Noise standard deviations of unit vector measurements
sigma_r = 1e-1;

% % Unit conversion
% P_a = P_a*(pi/180)^2; % rad^2
% P_g = P_g*(pi/(180*3600))^2; % (rad/sec)^2
% b0  = b0*pi/(180*3600); % rad/sec

% Variables
w_true = [w1_true; w2_true; w3_true];
P = [P_a, zeros(3); zeros(3), P_g];
dx_hat = zeros(6, length(time));
q_true = zeros(4, length(time));
q_hat  = zeros(4, length(time));
b_hat  = zeros(3, length(time));
w_hat  = zeros(3, length(time));
P_diag = zeros(6, length(time));
q_true(:,1) = q0;
q_hat(:,1)  = q0;
b_hat(:,1)  = b0;
dx_hat(:,1) = dx0;
sigma_g = diag([sigma_gx, sigma_gy, sigma_gz]);
sigma_b = diag([sigma_bx, sigma_by, sigma_bz]);

% Orientation in quaternions as the body rotates
for k = 1:length(time)-1
    q_true(:,k+1) = propQuaternion(q_true(:,k), w_true(:,k), dt); 
end

% Simulate gyro and unit vector measurements
[w_gyro, b_true] = simulateGyro(w_true, sigma_g, sigma_b, dt);
[v_b, v_r] = generateVectorPairs(q_true, sigma_r);

for k = 1:length(time)-1
    % Diagonal elems of P for 3-sigma bound
    P_diag(:,k) = diag(P);
    
    % Gyro-bias correction
    w_hat(:,k) = w_gyro(:,k) - b_hat(:,k); 
    
    % MEKF
    [q_hat(:,k+1), b_hat(:,k+1), dx_hat(:,k+1), P] = ... 
        mekf_murrell(q_hat(:,k), b_hat(:,k), w_hat(:,k), P, ...
        v_b(:,k), v_r(:,k), sigma_r, sigma_g, sigma_b, dt);
end 

% Sigma bound analysis
[sgma_bnd, q_err] = sigmaBoundErrors(q_true, q_hat, P_diag);
plotOutput(time, q_true, q_hat, dx_hat, sgma_bnd, ...
            w_hat, w_gyro, w_true, b_true, b_hat)


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
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

% Discrete-time gyro measurement simulation
function [w_gyro, bias] = simulateGyro(w_true, sigma_g, sigma_b, dt)
%
% Inputs:
%   sigma_g = Gyro noise standard deviation matrix [3x3]
%   sigma_b = Gyro-bias noise standard deviation matrix [3x3]
[~,n] = size(w_true);
sg1 = sigma_g(1,1); sg2 = sigma_g(2,2); sg3 = sigma_g(3,3);
sb1 = sigma_b(1,1); sb2 = sigma_b(2,2); sb3 = sigma_b(3,3);

% Gyro noise simulation
gn1 = sqrt(sg1^2/dt + sb1^2*dt/12)*randn(1,n); 
gn2 = sqrt(sg2^2/dt + sb2^2*dt/12)*randn(1,n); 
gn3 = sqrt(sg3^2/dt + sb3^2*dt/12)*randn(1,n);

% Gyro-bias noise simulation
bn1 = sb1/sqrt(dt)*randn(1,n);
bn2 = sb2/sqrt(dt)*randn(1,n);
bn3 = sb3/sqrt(dt)*randn(1,n);
gyro_noise = [gn1; gn2; gn3];

% Gyro-bias simulation
t = 0.1*pi/180/3600/dt;
[A, B, C, D] = tf2ss(dt*[1 1], 2*[1 -1]);
b1 = dlsim(A, B, C, D, bn1, t)';
b2 = dlsim(A, B, C, D, bn2, t)';
b3 = dlsim(A, B, C, D, bn3, t)';
bias = [b1; b2; b3];

% Gyro measurement
w_gyro = w_true + gyro_noise + bias;
end

% Generate measurement and reference vector pairs with 
% white noise of given standard-deviation.
function [v_b, v_r] = generateVectorPairs(q, sigma)
[~,m] = size(q);
v_r = rand(3,m);
v_b = zeros(3,m);
for i_iters = 1:m
    v_r(:,i_iters) = v_r(:,i_iters)/norm(v_r(:,i_iters));
    v_b(:,i_iters) = quaternion2A(q(:,i_iters))*v_r(:,i_iters);
end
% Zero-mean noise with sigma^2 variance 
noise = sigma.*randn(3,m);
v_b = v_b + noise;
end

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

function [sigma_bound, q_err] = sigmaBoundErrors(q_true, q_hat, P_diag)
% 3-sigma bound of states from diagonal covariances of the state error
% covariance matrix (Pg. 462)
sigma_bound = P_diag.^(0.5)*3*180/pi; % rad to deg

% Compute error quaternion using q_true and q_hat
% Note: 
%   Eq.(7.20) and Eq.(7.20) implies that the dx(1:3) = d_alpha are the 
%   components of roll, pitch & yaw error angles for any rotation sequence
q_hat(:,1:3) = -q_hat(:,1:3);
q_err(:,1) =    q_true(:,4).*q_hat(:,1) + q_true(:,3).*q_hat(:,2) ...
              - q_true(:,2).*q_hat(:,3) + q_true(:,1).*q_hat(:,4);
q_err(:,2) =  - q_true(:,3).*q_hat(:,1) + q_true(:,4).*q_hat(:,2) ... 
              + q_true(:,1).*q_hat(:,3) + q_true(:,2).*q_hat(:,4);
q_err(:,3) =    q_true(:,2).*q_hat(:,1) - q_true(:,1).*q_hat(:,2) ... 
              + q_true(:,4).*q_hat(:,3) + q_true(:,3).*q_hat(:,4);
q_err(:,4) =  - q_true(:,1).*q_hat(:,1) - q_true(:,2).*q_hat(:,2) ...
              - q_true(:,3).*q_hat(:,3) + q_true(:,4).*q_hat(:,4);
end

function plotOutput(time, q_true, q_hat, dx_hat, sgma_bnd, ...
            w_hat, w_gyro, w_true, b_true, b_hat)
% Quaternions
figure;
subplot(4,1,1); hold on;
plot(time,q_true(4,:));
plot(time,q_hat(4,:));
title('Quaternion attitude estimation');
legend('Ground truth','Estimate');
ylabel('q_0');

subplot(4,1,2);  hold on;
plot(time,q_true(1,:));
plot(time,q_hat(1,:));
ylabel('q_0');
legend('Ground truth','Estimate');

subplot(4,1,3); hold on;
plot(time,q_true(2,:));
plot(time,q_hat(2,:));
ylabel('q_1');
legend('Ground truth','Estimate');

subplot(4,1,4); hold on;
plot(time,q_true(3,:));
plot(time,q_hat(3,:));
ylabel('q_3'); xlabel('Time (sec)');
legend('Ground truth','Estimate');

% Three sigma bound
figure;
subplot(3,1,1); hold on;
plot(time,dx_hat(1,:)*3*180/pi);
plot(time,sgma_bnd(1,:),'color','k');
plot(time,-sgma_bnd(1,:),'color','k');
ylabel('Roll');
title('Attitude error (deg) and 3\sigma bound');

subplot(3,1,2); hold on;
plot(time,dx_hat(2,:)*3*180/pi);
plot(time,sgma_bnd(2,:),'color','k');
plot(time,-sgma_bnd(2,:),'color','k');
ylabel('Pitch');

subplot(3,1,3); hold on;
plot(time,dx_hat(3,:)*3*180/pi);
plot(time,sgma_bnd(3,:),'color','k');
plot(time,-sgma_bnd(3,:),'color','k');
ylabel('Yaw'); xlabel('Time (sec)');

% Angular velocity
figure;
subplot(3,1,1); hold on;
plot(time,w_hat(1,:),'.'); 
plot(time,w_gyro(1,:),'.');
plot(time,w_true(1,:),'LineWidth',1.5);
legend('Estimate','Gyro reading','Ground truth');
title('Angular velocity correction');
ylabel('\omega_1');

subplot(3,1,2); hold on;
plot(time,w_hat(2,:),'.'); 
plot(time,w_gyro(2,:),'.');
plot(time,w_true(2,:),'LineWidth',1.5);
legend('Estimate','Gyro reading','Ground truth');
ylabel('\omega_2');

subplot(3,1,3); hold on;
plot(time,w_hat(3,:),'.'); 
plot(time,w_gyro(3,:),'.');
plot(time,w_true(3,:),'LineWidth',1.5);
legend('Estimate','Gyro reading','Ground truth');
ylabel('\omega_3'); xlabel('Time (sec)');

% Bias
figure;
subplot(3,1,1); hold on;
plot(time,b_true(1,:));
plot(time,b_hat(1,:));
legend('True bias','Bias estimate');
title('Bias estimation (deg/hr)');
ylabel('x');

subplot(3,1,2); hold on;
plot(time,b_true(2,:));
plot(time,b_hat(2,:));
legend('True bias','Bias estimate');
ylabel('y');

subplot(3,1,3); hold on;
plot(time,b_true(3,:));
plot(time,b_hat(3,:));
legend('True bias','Bias estimate');
ylabel('z'); xlabel('Time (sec)');
end
