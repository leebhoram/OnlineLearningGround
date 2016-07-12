% by Bhoram Lee 
%
% Reference
% Titterton & Weston, "Strapdown Inertial Navigation Technology", 2nd. ed.,
% 2004.
function [ eulr ] = dcm2eulr( dcm )

N = size(dcm,3);
eulr = zeros(3,N);
for k=1:N
    eulr(2,k) = atan2(-dcm(1,3,k),sqrt(dcm(3,2,k)^2 + dcm(3,3,k)^2));
    
    if abs(eulr(2,k) - pi/2) < 10e-4 % near singularity
        eulr(1,k) = 0;
        eulr(3,k) = 0;   
    else
        eulr(1,k) = atan2(dcm(2,3,k),dcm(3,3,k));
        eulr(3,k) = atan2(dcm(1,2,k),dcm(1,1,k));    
    end
end

end

