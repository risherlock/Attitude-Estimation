function [q_hat, P_hat, beta_hat] = soar(q, P, Q, beta, v_i, v_b, w, dt)
% Inputs:
%   q_   = A-priori attitude estimate
%   P_   = A-priori covariance matrix
%   v_i  = Vector in inertial frame
%   v_b  = Vector sensor readings
%   w    = Gyroscope readings
%   dt   = Sample time
%   Q    = Process noise covariance [6x6] 
%
% Reference:
%   Sequential Optimal Attitude Recursion Filter
%   John A. Christian, E. Glenn Lightsey (2010)
%
% Rishav (2022/03/15)


% The algorithms proceeds in 4 steps:
%   A. Compute a-priori Davenpport matrix
%   B. Compute measurement Davenport matrix
%   C. Update state vector and covariance matrix
%   D. Propagate state estimate and covariance matrix
%
% Notes:
%   Propagated state estimate and covariance matrix of last time step is input 
%   to current time step. So, it is reasonable to perform this step at first.

% STEP D %
% Propagate quaternion and state covariance
q_ = propagate_quaternion(q, w, dt);
P_ = propagate_state_covariance(P, Q, w, dt);

% STEP A %
% Compute a-priori attitude profile matrix, B_
I = eye(3);
inv_Ptt = inv(P_(1:3, 1:3));
T = quaternion_to_dcm(q_);
B_ = (0.5*trace(inv_Ptt)*I - inv_Ptt)*T;

% Compute a-priori Davenport matrix, K_
S = B_ + B_';
mu = trace(B_);
z = [B_(2,3)-B_(3,2); B_(3,1)-B_(1,3); B_(1,2)-B_(2,1)];
K_ =[S - mu*I, z; z', mu];

% STEP B %
% Compute measurement attitude profile matrix, B_m
B_m = v_b*v_i';

% Compute measurement Davenport matrix, K_m
S = B_m + B_m';
mu = trace(B_m);
z = [B_m(2,3)-B_m(3,2); B_m(3,1)-B_m(1,3); B_m(1,2)-B_m(2,1)];
K_m =[S - mu*I, z; z', mu];

% STEP C %
% Compute SOAR attitude profile and Davenport matrix
B_hat = B_m + B_;
K_hat = K_m + K_;

% Update  quaternion attitude
q_hat = [1, 0, 0, 0]';

% Compute Ficher information matrix
F_ = inv(P_);
F_bt = F_(4:6, 1:3);
F_tb = F_(1:3, 4:6);
F_bb = F_(4:6, 4:6);

% Update non attiude parameters
Psi = [q_(4),q_(3),-q_(2);-q_(3),q_(4),q_(1);q_(2),-q_(1),q_(4);-q_(1),-q_(2),-q_(3)]';
dbeta = -2 * inv(F_bb) * F_bt * Psi * q_hat;
beta_hat = beta + dbeta;

% Update Ficher information matrix
Ftb = F_tb;
Fbt = Ftb';
Fbb = F_bb;
Ftt = trace(T * B_hat') * I - T * B_hat' + F_tb * inv(F_bb) * F_bt;
F_hat = [Ftt , Ftb; Fbt, Fbb];

% Update covariance matrix
P_hat  = inv(F_hat);
end
 