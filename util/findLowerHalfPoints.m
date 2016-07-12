function [ Lidx, Hidx, u_, v_] = findLowerHalfPoints(uv)

global cv
global cu

Hidx = find(uv(:,2)<=cv);
Lidx = find(uv(:,2)>cv);
v_ = uv(Lidx,2) - cv;
u_ = uv(Lidx,1) - cu;
end

