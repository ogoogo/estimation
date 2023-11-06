function [r] = christian_robinson(T_PB,s,Ap)
%CHRISTIAN_ROBINSON Performs the Christian-Robinson algorithm to calculate
%                   the displacement of a sensor to a celestial medium
%                   center
% 
% r = christain_robinson(q,s,Ap)
%
% q = rotation quaternion from body to planetary frame (q0,q1,q2,q3)
% s = horizon limb unit vector measurements (sx,sy,sz)
% Ap = ellipsoidal body matrix with diag(1/a2,1/b2,1/c2)
% sigma = noise of each horizon measurement
% d = pixel pitch in terms of pixels per radian
%
% r = estimated position relative to celestial body center in planet frame
%     (x,y,z)
%
% Author: Joshua Critchley-Marrows (USYD), 21-Feb-2022
% Modified: Joshua Critchley-Marrows (USYD), 2-Aug-2022

m = length(s);

% Specify the transformation matrix
% T_PB = [q(1)^2+q(2)^2-q(3)^2-q(4)^2 2*(q(2)*q(3)+q(1)*q(4))     2*(q(2)*q(4)-q(1)*q(3));
%         2*(q(2)*q(3)-q(1)*q(4))     q(1)^2-q(2)^2+q(3)^2-q(4)^2 2*(q(3)*q(4)+q(1)*q(2));
%         2*(q(2)*q(4)+q(1)*q(3))     2*(q(3)*q(4)-q(1)*q(2))     q(1)^2-q(2)^2-q(3)^2+q(4)^2];

% Calculate Christian-Robinson parameters
D = sqrt(Ap);
H = zeros(m,3);
for i = 1:m
    H(i,:) = (D*T_PB*s(:,i))';
    H(i,:) = H(i,:)/norm(H(i,:));
end
Dinv = zeros(3); 
Dinv(1,1) = 1/D(1,1); Dinv(2,2) = 1/D(2,2); Dinv(3,3) = 1/D(3,3);

% Perform TLS to find n
[~,~,V] = svd([H ones(m,1)],0);
n = -V(1:3,4)/V(4,4);
if H(1,1:3)*n < 0, n = -n; end
% n = (H'*H)\H'*ones(m,1);

% Calculate position in planetary frame
r = -1/sqrt(n'*n-1)*Dinv*n;

% Compute the variation and covariance

% Rs = (sigma/d)^2*[1,0,0;0,1,0;0,0,0];
% Ry = zeros(m);
% for i = 1:m
%     sbar = D*T_PB*s(:,i);
%     J = n'/norm(sbar)*(eye(3) - sbar*sbar'/norm(sbar)^2);
%     Ry(i,i) = J*D*Rs*D'*J';
% end
% Pn = pinv(H'*pinv(Ry)*H);
% F = -Dinv/sqrt(n'*n-1)*(eye(3) - n*n'/(n'*n-1));
% P = F*Pn*F';
% d = sqrt(norm(eig(P)));
% 
end