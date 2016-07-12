function [bestModel, inliers, distSqr, ransac_info]  = computeScale_3param(data,weights)
% data := [XinLowerHalf/f/h YinLowerHalf/f, Zinv_inLowerHalf, triMatch(groundIdx,7:8)]

%bestModel = [0 0 0];   
%inliers = []; % indices of ground inliers

% scale data to make them roughly range between [-1 1]
x = 2*data(:,1);   
y = 2*data(:,2);   
med3 = 3*median(data(:,3));
d = data(:,3)/med3;   

N = size(data,1);
A = [x.*d y.*d ones(N,1)];
b = d.*d;

% bestModel = pinv(A)*b;

[bestModel, inliers, distSqr, ransac_info]= ransac_2ndPoly_3param(A,b,x,y,d, weights); 
bestModel = bestModel(2)*med3*2;
end % end of function

 
function [final_model, final_inliers, distSqr, ransac_info] = ransac_2ndPoly_3param(A,b,x,y,d,weights)
   
global Ransac_param

k = 0;              % count of saved hypothesis
ite = 0;            % count of iteration
model = [];
final_inliers = 0;
final_model = [0 0 0];
N = length(b);
last_num_cons = 0;
% err = NaN;
ransac_info = struct('last_ite', 0, 'inlierRatio', 0,'maxllhModel',[0 0 0]);

if ~isempty(weights)
    AccWeights = weights; % accumulated weights
    for k=2:N
        AccWeights(k) = AccWeights(k-1) + weights(k);
    end
    AccWeights = AccWeights/AccWeights(end);
end


% required_num_data = round(Ransac_param.req_numData*N);
data_thre =  Ransac_param.sig_dist^2;
%figure(3), hold off; plot3(A(:,1), A(:,2), b, '.'); hold on;
        

if N > 3
    while ite < Ransac_param.req_numIter 
        if ~isempty(weights)
            idx = ones(3,1);
            % weighted sampling      
            while length(unique(idx))<3
                rndnum = rand([3 1]);
                for s = 1:3            
                    idx(s) = findCeilingIndex(AccWeights,rndnum(s));
                   % disp(strcat('Rand#:',num2str(rndnum(s)),', index:',int2str(idx(s))));
                end
            end
        else
            idx = randi(N,[3 1]);
        end
        
        A_ = A(idx,:);
        b_ = b(idx);
        %figure(3), plot3(A_(:,1), A_(:,2), b_, 'ro','MarkerFaceColor','r'); hold on;
        
        if cond(A_) < 1e3
            m_ = A_\b_;
            allDist = (A*m_ - b).^2;            
            inlierIdx = find(allDist<data_thre);
                  num_cons = length(inlierIdx);
            
            if num_cons > 3 
                A_ = A(inlierIdx,:);
                b_ = b(inlierIdx);
                model = A_\b_;
                                
                if (k == 0) || (num_cons > last_num_cons)
                     k = k + 1;
                     last_num_cons = num_cons;
                     final_inliers = inlierIdx;
                     final_model = model;
                     
                     ransac_info.last_ite = ite;
                     ransac_info.inlierRatio = num_cons/N;
                end                             
            end            
        end
        
        ite = ite + 1;
    end
end

ransac_info.maxllhModel = final_model;  
distSqr = (A*final_model - b).^2; 

% figure(13), 
% plot3(x, y, d, 'k.'); hold on;
% plot3(x(final_inliers), y(final_inliers), d(final_inliers), 'ro','MarkerFaceColor','k');
% 
% figure(14), 
% plot3(A(:,1), A(:,2), b, 'k.'); hold on;
% plot3(A(final_inliers,1), A(final_inliers,2), b(final_inliers), 'ro','MarkerFaceColor','k');
     
end


function idx = findCeilingIndex(orderedData,target)

diff = orderedData - target;
idx = min(find(diff>0));
end

