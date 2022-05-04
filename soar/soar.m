function [q_hat] = soar(q, w, P, v_i, v_b)
% Implementation of SOAR filter
% 
% References:
%   [1] Sequential Optimal Attitude Recursion Filter
%       John A. Christian, E. Glenn Lightsey (2010)
%
% Rishav (2022/03/15)

% [Step 1] Compute a-priori Davenport matrix
F = inv(P); % Eqn(61)
F_tt = F(1:3,1:3); 
T = quaternion_to_dcm(q);
B = (0.5*trace(F_tt)*eye(3) - F_tt) * T; % Eqn(33) 

% Construct Davenport matrix, Eqn(12) and Eqn(13)
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
K_minus = [B + B'- eye(3)*trace(B), Z; Z', trace(B)];

% [Step 2] Compute measurement Davenport matrix
B_m = (v_b.*repmat(w,[1 3])')*v_i';
Z_m = [B_m(2,3)-B_m(3,2); B_m(3,1)-B_m(1,3); B_m(1,2)-B_m(2,1)];
K_m = [B_m + B_m'- eye(3)*trace(B_m), Z_m; Z_m', trace(B_m)];

% [Step 3] Update state vector & covariance Matrix
K_plus = K_m + K_minus;
Psi = 0;

q_plus = esoq();


% [Step 4] Update state estimate and covariance to next measurement time
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
