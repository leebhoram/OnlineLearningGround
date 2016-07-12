% by Bhoram Lee 
%
% Reference
% Titterton & Weston, "Strapdown Inertial Navigation Technology", 2nd. ed.,
% 2004.
function [ dcm ] = quat2dcm(quat)

N = size(quat,2);
dcm = zeros(3,3,N);
for k=1:N
   q = quat(:,k);
   dcm(1,1,k) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4);
   dcm(1,2,k) = 2*(q(2)*q(3) + q(1)*q(4));
   dcm(1,3,k) = 2*(q(2)*q(4) - q(1)*q(3));
   
   dcm(2,1,k) = 2*(q(2)*q(3) - q(1)*q(4));
   dcm(2,2,k) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4);
   dcm(2,3,k) = 2*(q(1)*q(2) + q(3)*q(4));
   
   dcm(3,1,k) = 2*(q(2)*q(4) + q(1)*q(3));
   dcm(3,2,k) = 2*(q(3)*q(4) - q(1)*q(2));
   dcm(3,3,k) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4);     
end

end

