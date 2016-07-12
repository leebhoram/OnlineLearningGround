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
function mseq = matchFeature(pts1,pts2,fvec1,fvec2,mseq,V,alpha)

global f cu cv dt width height 

% ==== Update Inverse Depth
ZinvQ = 0.02; % inverse depth noise variance 
zPrevIdx = find(mseq.invDepthFlag > 0);

if zPrevIdx > 0
    % correct the update equation
    for tt = 1:length(zPrevIdx)
        d0 = mseq.invDepth(zPrevIdx(tt));
        uk = double((mseq.m1(zPrevIdx(tt),3)-cu))/f;
        z1 = double((1/d0*(1+uk*alpha*dt) - V*dt));
        if z1 > 0
            d1 = 1/z1;
            mseq.invDepthPred(zPrevIdx(tt)) = d1;    
            J = (d1/d0)^2*(1+uk*alpha*dt);
            mseq.invDepthVar(zPrevIdx(tt)) = J*J*mseq.invDepthVar(zPrevIdx(tt)) + ZinvQ;
        end
    end
end

% ==== match
if ~isempty(mseq.invDepthPred)
    invDepth1 = zeros(length(pts2),1);
    invDepth1(mseq.m1(:,6)) = mseq.invDepth;
else
    invDepth1 = zeros(length(pts2),1);
end 

% predict where to find each corresponce
Rsearch = 5;
param_temp = [f cv cu dt width height Rsearch alpha V];
[ pts1_pred, ~,~] = mexPredictBackFlow( double(pts2(:,1:2)'), [], param_temp);
% and match

SAD_thre = 500;
binSize = 40;
nxbin = ceil(width/binSize);
nybin = ceil(height/binSize);
FMparam = double([2 alpha width height  dt f cu cv Rsearch SAD_thre binSize nxbin nybin]);    
[matches] = mexMatchFeatureB_Geiger(pts2, pts1, fvec2, fvec1,int16(pts1_pred) ,FMparam);
% write in mseq
%                   (u,v)image2     (u,v)image1    indx(2nd)     indx(1st)
mseq.m2 = [matches(:,5:6) matches(:,3:4) matches(:,2) matches(:,1) ];    

% init mseq if needed
if isempty(mseq.invDepth)
    n2 = length(mseq.m2); 
    mseq.age = zeros(n2,1); 
    mseq.gndLLH = zeros(n2,1); 
    mseq.invDepth = zeros(n2,1); 
    mseq.invDepthPrev = zeros(n2,1); 
    mseq.invDepthFlag = zeros(n2,1); 
    mseq.invDepthPred = zeros(n2,1); 
    mseq.invDepthEst = zeros(n2,1); 
    mseq.invDepthVar = 2*ones(n2,1);
end

end