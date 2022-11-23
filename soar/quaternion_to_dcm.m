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
