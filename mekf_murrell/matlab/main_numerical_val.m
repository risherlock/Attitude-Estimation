%%% Numerical validation with C++ implementation.
% (Check the file ../cpp/main_mekf.cpp)
% 2020/1/18
clc
clear
close all

q = [0,0,0,1]'; % Initial quaternion
b = [1.324,2.4532,3.4532]'; % Initial bias
P = eye(6); % Initial state noise covariance
dt = 0.134; % Time step

% Noise standard deviation
sigma_r = [1.345, 2.654];
sigma_g = diag([4.453, 2.534, 3.234]);
sigma_b = diag([2.435, 1.2345, 4.55]);

% Vector measurements and reference vectors
mb = [1.32,2.65,3.63; 4.6534,5.634,6.634]';
mr = [7.53,8.23,9.97; 1.42,2.432,3.546]';

% Gyro reading 
w = [0.223,2.34,3.45]';

% Run EKF
[q_hat, b_hat, dx_hat, P_k] = ...
mekf_murrell(q, b, w - b, P, mb, mr, sigma_r, sigma_g, sigma_b, dt);

% Display
fprintf("Estimated quaternion:\n"); display(q_hat');
fprintf("Estimated bias:\n"); display(b_hat');
