% Barebones implementation of SOAR algorithm

% SETUP: This part of code runs only once %
% Initial conditions %
dt = 1;             % Sample time, seconds
q = [1, 0, 0, 0]';  % Initial quaternion
P = eye(6);         % Initial state covariance matrix
Q = eye(6);         % Process noise covariance matrix
beta = [0, 0, 0]';  % Initial gyro drift 

% LOOP: This part of code runs for each sensor readings %
% Sensor readings %
m_b = [1,2,3]'; m_i = [1,2,3]'; % Magnetometer
s_b = [1,2,3]'; s_i = [1,2,3]'; % Sun sensor
w = [1, 2, 3];                  % Gyroscope

% Vectors in body and inertial frame
v_b = [m_b, s_b];
v_i = [m_i, s_i];

% Gyro bias correction
w = w - beta;

% Run SOAR in loop
[q_hat, P_hat, beta_hat] = soar(q, P, Q, beta, v_i, v_b, w, dt);
q = q_hat;
P = P_hat;
beta = beta_hat;

printf("Quaternion estimate:\r\n");
display(q_hat);
printf("State covariance matrix:\r\n");
display(P_hat);
printf("Gyro drift:\r\n");
display(beta_hat);
