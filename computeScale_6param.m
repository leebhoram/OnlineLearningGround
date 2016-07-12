% Weighted Regularized Least Squares
function [bestModel, inliers]  = computeScale_6param_v0_demo(data, weights)

y = 2*data(:,2);   
x = 2*data(:,1); 
med3 = 3*median(data(:,3));
b = data(:,3)/med3;  % inverse depth of arbitrary scale

A = [ y.*b  x  b  x.*x  x.*b ones(length(x),1)];
B = b.*b;
[bestModel, inliers ] = nonflat_fit_weightedRansac(A,B, weights); 

bestModel = bestModel(1)*med3*2;

end % end of function


function [final_model, inliersFinal] = nonflat_fit_weightedRansac(A,b,weights)

global Ransac_param6;

N = length(b);
AccWeights = weights; % accumulated weights

for k=2:N
    AccWeights(k) = AccWeights(k-1) + weights(k);
end

AccWeights = AccWeights/AccWeights(end);
nParam = size(A,2);

k = 0;
model = zeros(1,nParam);
inliersFinal = [];

% criterion
required_itr = Ransac_param6.req_numIter;
data_thre = Ransac_param6.sig_dist.^2;
final_num = 0;
itr = 1;
Ite = 0;
k = 0;


samples = [];
Ns = N;
rng('shuffle');  
while itr < required_itr && N > 1
    
     if ~isempty(weights)
            idx = zeros(nParam,1);
            % weighted sampling      
            
            for s = 1:nParam               
                if s == 1 
                    rndnum = rand(1);
                    idx(s) = findCeilingIndex(AccWeights,rndnum);
                else
                    idx(s) = idx(s-1);
                    iit = 0;
                    while ~isempty(find(idx(1:(s-1)) == idx(s))) && iit < 30
                         rndnum = rand(1);
                        idx(s) = findCeilingIndex(AccWeights,rndnum);
                        iit = iit + 1;
                    end
                end
            end
    else
        idx = randi(N,[nParam 1]);
    end
        
    samples = [samples; idx];
    A_ = A(idx,:);
      
    if cond(A_) < 1e5               
      
        A__ = A_'*diag(weights(idx))*A_ + 0.1*diag([0 1 1 1 1 1]);
        bb = A__\(A_'*diag(weights(idx))*b(idx));
      
        % consensus check
        allDist = (A*bb - b)/sqrt(sum(bb(1:(end-1)).^2 + 1)); 
        allSqrDist = allDist.*allDist;
        err = sqrt(mean(allSqrDist));
        InliersInd = find(allSqrDist < data_thre);
        num_cons = length(InliersInd);
        
        if  (num_cons >= 2) % we may have found a good model

            % best? 
            if k == 0 || (num_cons >  final_num(k))
                k = k + 1;
                final_num(k) = num_cons;
                model(k,:) = bb;
                inliersFinal = InliersInd;  
                Dist = allDist;           
                Ite(k) = itr;
            end
        end
    end
    
    itr = itr + 1;
end

if k == 0
    disp('!!! Fail to find fit'); 
    inliersFinal = 0;    
    [ bb] = iteReWeightedLS2( A,b,[],1);
    model = bb;
    final_model = model;
    Dist = A*bb - b; 
    final_num = 0;
    
else
   % final refinement
    A_ = A(inliersFinal,:);
    A__ = A_'*diag(weights(inliersFinal))*A_ + 0.1*diag([0 1 1  1 1 1]);
    bb = A__\(A_'*diag(weights(inliersFinal))*b(inliersFinal));
    
    model(end,:) = bb;
    final_model = model(end,:);
end


end


function idx = findCeilingIndex(orderedData,target)

diff = orderedData - target;
idx = min(find(diff>=0));
end

