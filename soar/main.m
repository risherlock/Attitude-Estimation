% Test of SOAR algorithm
% 
% References:
%   Sequential Optimal Attitude Recursion Filter
%   John A. Christian, E. Glenn Lightsey (2010)
%
% Rishav (2022/03/15)

clc
clear
close all

% Simulation params
dt = 0.5; % Step size, s
stop_time = 100; % End time, s
time = 0:dt:stop_time;

% Angular velocities (ground truth)
wx_true = 0.1 + 0.1*sin(time+1); % Rotation about x-axis, rad/s
wy_true = 0.2 + 0.3*sin(time+2); % Rotation about y-axis, rad/s
wz_true = 0.3 + 0.9*sin(time+3); % Rotation about z-axis, rad/s

% Quaternions: [q0, q1, q2, q3] -> [cos(psi/2), sin(psi/2) * v]
q0 = [1, 0, 0, 0]'; % Initial quaternion

% Sensor reading in inertial frame
m_i = [0.6651, 0.7395, 0.1037]'; % Unit vector compass reading
s_i = [0.8190, 0.5671, 0.0874]'; % Unit vector sun vector reading

% Intermediate variables
m_b = zeros(3, length(time)); % Compass reading in body frame
s_b = zeros(3, length(time)); % Sun vector reading in body frame
q_true = zeros(4, length(time));
w_true = [wx_true; wy_true; wz_true];
q_true(:,1) = q0;

% Orientation in quaternions as the body rotates
for k = 1:length(time)
    if(k ~= length(time)) 
        q_true(:,k+1) = propagate_quaternion(q_true(:,k), w_true(:,k), dt);
    end
    Q = quaternion_to_dcm(q_true(:,k)); % Rotation matrix
    m_b(:,k) = Q * m_i;
    s_b(:,k) = Q * s_i;
end

figure;
plot(time,q_true(1,:)); hold on;
plot(time,q_true(2,:));
plot(time,q_true(3,:));
plot(time,q_true(4,:)); grid on;
legend('q_0','q_1','q_2','q_3');
title('Quaternion profile');

figure;
plot(time,w_true(1,:)); hold on
plot(time,w_true(2,:));
plot(time,w_true(3,:)); grid on;
legend('w_x','w_y','w_z');
title('Angular velocity profile');

figure;
plot(time,m_b(1,:)); hold on
plot(time,m_b(2,:));
plot(time,m_b(3,:)); 
plot(time,(m_b(1,:).^2 + m_b(2,:).^2 + m_b(3,:).^2).^0.5); grid on;
legend('m_x','m_y','m_z','norm');
title('Magnetic field profile in body frame');

figure;
plot(time,s_b(1,:)); hold on
plot(time,s_b(2,:));
plot(time,s_b(3,:));
plot(time,(m_b(1,:).^2 + m_b(2,:).^2 + m_b(3,:).^2).^0.5); grid on;
legend('s_x','s_y','s_z','norm');
title('Sun vector profile in body frame');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Discrete-time quaternion propagation
function [q_next] = propagate_quaternion(q, w, dt)
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
    q_next = Omega*q; % Eq.(7.39)
else
    q_next = q;
end
end

% Quaternion to rotation matrix
function [Q] = quaternion_to_dcm(q)
% Input:
%   q = [q0,q1,q2,q3]
%     = [cos(psi/2), e1 * sin(psi/2), e2 * sin(psi/2), e3 * sin(psi/2)]
%   
%   where, 
%       e1,e2,e3 = 3D orthogonal basis vectors and
%       psi = Angle of the axis-angle representation of quaternion
%
% Output:
%   Q = Rotation matrix corresponding to the input quaternion
%
% Reference:
%   [1] Schaub, Junkins - Analytical Mechanics of Space Systems (4th ed.)

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

% Eqn(3.98) [1]
Q = zeros(3);
Q(1,1) = q0^2 + q1^2 - q2^2 - q3^2;
Q(2,2) = q0^2 - q1^2 + q2^2 - q3^2;
Q(3,3) = q0^2 - q1^2 - q2^2 + q3^2;
Q(1,2) = 2*(q1*q2 + q0*q3);
Q(1,3) = 2*(q1*q3 - q0*q2);
Q(2,1) = 2*(q1*q2 - q0*q3);
Q(2,3) = 2*(q2*q3 + q0*q1);
Q(3,1) = 2*(q1*q3 + q0*q2);
Q(3,2) = 2*(q2*q3 - q0*q1);
end
