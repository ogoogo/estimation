function [rp,q] = nlls_horizon(rp,q,s,Ap)
%NLLS_HORIZON Computes an NLLS solution for attitude and position based on 
%             a horizon estimate
% 
% [rp,q] = nlls_horizon(rp,q,s,Ap)
%
% rp = guessed position relative to celestial body center in planet frame
%     (x,y,z)
% q = rotation quaternion from body to planetary frame (q0,q1,q2,q3)
% s = horizon limb unit vector measurements (sx,sy,sz)
% Ap = ellipsoidal body matrix with diag(1/a2,1/b2,1/c2)
%
% rp = estimated position relative to celestial body center in planet frame
%     (x,y,z)
% q = estimated quaternion (body->inertial)
%
% Author: Joshua Critchley-Marrows (USYD), 24-Feb-2022
% Modified: Joshua Critchley-Marrows (USYD), 1-Mar-2022
%

m = length(s);
thresh = 10;
iter_max = 100;

% Compute model parameters
T_PB = [q(1)^2+q(2)^2-q(3)^2-q(4)^2 2*(q(2)*q(3)+q(1)*q(4))     2*(q(2)*q(4)-q(1)*q(3));
        2*(q(2)*q(3)-q(1)*q(4))     q(1)^2-q(2)^2+q(3)^2-q(4)^2 2*(q(3)*q(4)+q(1)*q(2));
        2*(q(2)*q(4)+q(1)*q(3))     2*(q(3)*q(4)-q(1)*q(2))     q(1)^2-q(2)^2-q(3)^2+q(4)^2];
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
dTdq0 = [2*q0,2*q3,-2*q2;-2*q3,2*q0,2*q1;2*q2,-2*q1,2*q0]';
dTdq1 = [2*q1,2*q2,2*q3;2*q2,-2*q1,2*q0;2*q3,-2*q0,-2*q1]';
dTdq2 = [-2*q2,2*q1,-2*q0;2*q1,2*q2,2*q3;2*q0,2*q3,-2*q2]';
dTdq3 = [-2*q3,2*q0,2*q1;-2*q0,-2*q3,2*q2;2*q1,2*q2,2*q3]';
Xi = [-q1,-q2,-q3;q0,-q3,q2;q3,q0,-q1;-q2,q1,q0];
    Mp = Ap*(rp*rp')*Ap - (rp'*Ap*rp-1)*Ap;
A = T_PB'*Ap*T_PB;
r = T_PB'*rp;
M = A*(r*r')*A - (r'*A*r - 1)*A;

% Set up NLLS parameters
delta = thresh^2*ones(6,1);
H = zeros(m,6);
h = zeros(m,1);
iter = 0;
while norm(delta) > thresh && iter < iter_max

    % Construct Jacobian and model vector
    for i = 1:m
        H(i,1:3) = 2*r'*A*s(:,i)*(A*s(:,i))' - 2*(s(:,i)'*A*s(:,i))*r'*A;
        dhdq(1) = (dTdq0'*s(:,i))'*Mp*T_PB*s(:,i) + ...
                    (T_PB*s(:,i))'*Mp*dTdq0'*s(:,i);
        dhdq(2) = (dTdq1'*s(:,i))'*Mp*T_PB*s(:,i) + ...
                    (T_PB*s(:,i))'*Mp*dTdq1'*s(:,i);
        dhdq(3) = (dTdq2'*s(:,i))'*Mp*T_PB*s(:,i) + ...
                    (T_PB*s(:,i))'*Mp*dTdq2'*s(:,i);
        dhdq(4) = (dTdq3'*s(:,i))'*Mp*T_PB*s(:,i) + ...
                    (T_PB*s(:,i))'*Mp*dTdq3'*s(:,i);
        H(i,4:6) = (0.5*dhdq*Xi);
        h(i) = s(:,i)'*M*s(:,i);
    end

    % Compute error estimate and update
    delta = -pinv(H'*H)*H'*h;
    r = r + delta(1:3);
    q = q + 0.5*Xi*delta(4:6);
    q = q/norm(q);

    % Recompute model parameters
    T_PB = [q(1)^2+q(2)^2-q(3)^2-q(4)^2 2*(q(2)*q(3)+q(1)*q(4))     2*(q(2)*q(4)-q(1)*q(3));
            2*(q(2)*q(3)-q(1)*q(4))     q(1)^2-q(2)^2+q(3)^2-q(4)^2 2*(q(3)*q(4)+q(1)*q(2));
            2*(q(2)*q(4)+q(1)*q(3))     2*(q(3)*q(4)-q(1)*q(2))     q(1)^2-q(2)^2-q(3)^2+q(4)^2];
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    dTdq0 = [2*q0,2*q3,-2*q2;-2*q3,2*q0,2*q1;2*q2,-2*q1,2*q0]';
    dTdq1 = [2*q1,2*q2,2*q3;2*q2,-2*q1,2*q0;2*q3,-2*q0,-2*q1]';
    dTdq2 = [-2*q2,2*q1,-2*q0;2*q1,2*q2,2*q3;2*q0,2*q3,-2*q2]';
    dTdq3 = [-2*q3,2*q0,2*q1;-2*q0,-2*q3,2*q2;2*q1,2*q2,2*q3]';
    Xi = [-q1,-q2,-q3;q0,-q3,q2;q3,q0,-q1;-q2,q1,q0];
    M = A*(r*r')*A - (r'*A*r - 1)*A;
    A = T_PB'*Ap*T_PB;
    rp = T_PB*r;
    Mp = Ap*(rp*rp')*Ap - (rp'*Ap*rp-1)*Ap;
    iter = iter + 1;
end
end