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
close all;
clear all;

imagefolder = '~/Datasets/kitti/odo_data/00/image_2/';
posefile = '~/Datasets/kitti/poses/00.txt';

addpath('./util');
addpath('./mex');
addpath('./display');

% ==== Load image list
[imagelist]= getImageListFromFolder(imagefolder);
N = length(imagelist);

% ==== Set Options
startIdx = 1;
endIdx = min(1000,length(imagelist));
mode = 3; % 6 : 6-param model for the ground, otherwise use 3-param model

% ==== Initialize
loadParams

% ====  the first image 
ImgOrig1 = impyramid(imread(imagelist{startIdx}),'reduce'); 
width = size(ImgOrig1,2);
height = size(ImgOrig1,1);
[IGrey1, ~ ,Pts1, Fvec1] = computeFeature(ImgOrig1); 

% ==== the Loop
for fn = (startIdx+1):(endIdx-1)
    % tic,
    ImgOrig2 = impyramid(imread(imagelist{fn}),'reduce');       
    
    % ==== feature computation and matching (based on Geiger's algorithm)
    [IGrey2, IndHSV,Pts2, Fvec2] = computeFeature(ImgOrig2); 
    MatchSequence = matchFeature(Pts1,Pts2,Fvec1,Fvec2,MatchSequence,Vhat,alphaHat);
    
    displayCurrentImage_fig1; 
    displayCurrentImage_fig2; 
    
    if isempty(MatchSequence.m1) 
         % ==== Move Forward     
        MatchSequence.m1 = MatchSequence.m2;   
        mseq.invDepthFlag = 0;
        Pts1 = Pts2;
        Fvec1 = Fvec2;
        continue;             
    else    
        % compute triplets (motion is computed based on central difference)
        triMatch = computeTriplet(MatchSequence);                
        displayTriMatches_fig1;                                

        if is_static(MatchSequence) == true
            V = 0; % force velocity to zero
            MatchSequence.invDepthFlag(:) = 0;  % no depth obvervation
            [VF, Vhat] = VF.zero();
            [wF, alphaHat] = wF.zero();
            disp('static!');
        else
            [MatchSequence, V, alpha, groundIdx_m1] = computeMotionAndScale(triMatch, MatchSequence, GndProb, mode); % NOTE: non-holonomic motion and fixed camera height assumption 
                    
            if isempty(groundIdx_m1) || length(groundIdx_m1) < 10
                [VF, Vhat] = VF.propagate();
                [wF, alphaHat] = wF.update(alpha);
            else                
                % motion filter update
                if ~isempty(Vhat) &&  abs(V-Vhat) > 1.7321
                    [VF, Vhat] = VF.update(V,((V-Vhat))^2); % update adaptively 
                else
                    [VF, Vhat] = VF.update(V);
                end
                [wF, alphaHat] = wF.update(alpha);                                                     

                % Update Inverse Depth
                MatchSequence = updateInvDepth(MatchSequence) ;                                                                    

                % Color learning
                [GndProb, MatchSequence, HSVbinsNG, HSVbinsG ] = learnGround(IndHSV,MatchSequence,groundIdx_m1,HSVbinsNG,HSVbinsG);

            end % end of if bVelocity_estimated == 1 or 0

        end % end of if static == 1 or 0

        disp(strcat('V_hat:', num2str(Vhat))); 
        yawHat = yawHat - alphaHat*dt;
        P = P + dt*[Vhat*sin(-yawHat) Vhat*cos(-yawHat)];
        Phat_mat = [Phat_mat; fn P]; % save trajectory 
        % toc
        
        disp('Press a key to continue.');
        pause; 
              
        % ==== Move Forward     
        MatchSequence.m1 = MatchSequence.m2;   
        updateMatchSequence;
        Pts1 = Pts2;
        Fvec1 = Fvec2;
        IGrey2 =[];
        Pts2 = [];
        Fvec2 = [];
        matches = [];
    end        

end

% trajectory plot 
figure(99), 
plot(ref_path(:,1), ref_path(:,3),'r-'); hold on;
plot(Phat_mat(:,2), Phat_mat(:,3),'.-'); axis equal; 
