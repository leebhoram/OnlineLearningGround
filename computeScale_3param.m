function [bestModel, inliers]  = computeScale_3param(data,weights)

% scale data to make them roughly range between [-1 1]
x = 2*data(:,1);   
y = 2*data(:,2);   
med3 = 3*median(data(:,3));
d = data(:,3)/med3;   

N = size(data,1);
A = [x.*d y.*d ones(N,1)];
b = d.*d;

[bestModel, inliers]= ransac_2ndPoly_3param(A,b, weights); 
bestModel = bestModel(2)*med3*2;
end % end of function

 
function [final_model, final_inliers] = ransac_2ndPoly_3param(A,b,weights)
   
global Ransac_param

k = 0;              % count of saved hypothesis
ite = 0;            % count of iteration
model = [];
final_inliers = 0;
final_model = [0 0 0];
N = length(b);
last_num_cons = 0;

if ~isempty(weights)
    AccWeights = weights; % accumulated weights
    for k=2:N
        AccWeights(k) = AccWeights(k-1) + weights(k);
    end
    AccWeights = AccWeights/AccWeights(end);
end

data_thre =  Ransac_param.sig_dist^2;
        

if N > 3
    while ite < Ransac_param.req_numIter 
        if ~isempty(weights)
            idx = ones(3,1);
            % weighted sampling      
            while length(unique(idx))<3
                rndnum = rand([3 1]);
                for s = 1:3            
                    idx(s) = findCeilingIndex(AccWeights,rndnum(s));
                end
            end
        else
            idx = randi(N,[3 1]);
        end
        
        A_ = A(idx,:);
        b_ = b(idx);
        
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
                end                             
            end            
        end
        
        ite = ite + 1;
    end
end

end


function idx = findCeilingIndex(orderedData,target)

diff = orderedData - target;
idx = min(find(diff>0));
end

