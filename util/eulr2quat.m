% by Bhoram Lee 
%
% Reference
% Titterton & Weston, "Strapdown Inertial Navigation Technology", 2nd. ed.,
% 2004.
function [ quat ] = eulr2quat( eulr )

N = size(eulr,2);
quat = zeros(4,N);
for k=1:N
    quat(1,k) = cos(eulr(1)/2)*cos(eulr(2)/2)*cos(eulr(3)/2) + sin(eulr(1)/2)*sin(eulr(2)/2)*sin(eulr(3)/2);
    quat(2,k) = sin(eulr(1)/2)*cos(eulr(2)/2)*cos(eulr(3)/2) - cos(eulr(1)/2)*sin(eulr(2)/2)*sin(eulr(3)/2);
    quat(3,k) = cos(eulr(1)/2)*sin(eulr(2)/2)*cos(eulr(3)/2) + sin(eulr(1)/2)*cos(eulr(2)/2)*sin(eulr(3)/2);
    quat(4,k) = cos(eulr(1)/2)*cos(eulr(2)/2)*sin(eulr(3)/2) - sin(eulr(1)/2)*sin(eulr(2)/2)*cos(eulr(3)/2);
    quat(:,k) = quat(:,k)/norm(quat(:,k));
end

end

