% Copyright (c) 2015 Bhoram Lee
% 
% This file is part of OLG.
% Authors: Bhoram Lee
% 
% OLG is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
function [mseq, V, alpha, groundIdx_m1] = computeMotionAndScale(tri, mseq, Gprob)

global h        % camera height from the ground
global f dt

groundIdx_m1 = [];
V = 0;
alpha = 0;

if size(tri,1) < 30 % too few observation    
    mseq.invDepthFlag(:) = 0;
    mseq.invDepth(:) = 0;
else
    % ==== Estimate inverse depth and yaw increment
    % NOTE: this function may be replaced with any motion estimation and 3D
    % reconstruction algorithm    
    [ ~, zIdx, Zinv, alpha] = computeIncMotion( double(tri), 5);    
    mseq.invDepthFlag(tri(zIdx,7)) = 1;  
    
    displayInvDepthOfPoints_fig1_demo;

    [lowerIdx, upperIdx, XinLowerHalf, YinLowerHalf] = findLowerHalfPoints(tri(:,3:4));
    [zLowerIdx, An, Bn]= intersect(lowerIdx,zIdx);
    Zinv_inLowerHalf = Zinv(Bn); 
    zYinLowerHalf = YinLowerHalf(An);
    zXinLowerHalf = XinLowerHalf(An);

    % upper half points unlikely to be GROUND point!
    mseq.gndLLH(tri(upperIdx,7)) = mseq.gndLLH(tri(upperIdx,7)) - 1;

    lowerHalf = [zYinLowerHalf/f/h, Zinv_inLowerHalf, zLowerIdx, (zXinLowerHalf)/f/h];   
    weights = ones(length(zLowerIdx),1);
    if ~isempty(Gprob)     
        % If ground probability is available, then use the weight (weighted ransac)
        zLowerm1 = tri(zLowerIdx,7);            
        for tt=1:length(zLowerm1)
            weights(tt) = Gprob(mseq.m1(zLowerm1(tt),4), mseq.m1(zLowerm1(tt),3))^2;
        end       
    end

    % compute the scale by fitting ground surface model (3 param) 
    [V, inliers, ~, ~] = computeScale_3param(lowerHalf(:,[4 1 2]),weights); 

     
    mseq.invDepth(tri(zIdx,7)) = Zinv/V; 
    
    if length(inliers) >= 10
        triGroundIndex = lowerHalf(inliers,3); % inlier_zLowerIdx of triMatch                
        groundIdx_m1 = tri(triGroundIndex,7);
        groundIdx_m2 = tri(triGroundIndex,8);
        mseq.gndLLH(groundIdx_m1) = mseq.gndLLH(groundIdx_m1) + 1;
    end
        
    display_v_d
    displayGroundPoints_fig1
    % NOTE: many outliers where d is close to zero.    
end


