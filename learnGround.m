function [GndProb, mseq, HSVbinsNG, HSVbinsG ]= learnGround(hsv,mseq,groundIdx,HSVbinsNG,HSVbinsG)

global width height
global Nquan

p_g = 0.2;      % prob(ground) (prior)
log_p_g = log(p_g);
p_ng = 1 - p_g; % prob(non-ground) (prior) 
log_p_ng = log(p_ng);
fadingRate = 0.7;
wndSampling = 1;
wndTesting = 1;

N_g = length(groundIdx);

accNGIndexm1 = intersect(find(mseq.gndLLH < 0),find(mseq.invDepthFlag == 1));

POI_NG = [mseq.m1(accNGIndexm1,3) mseq.m1(accNGIndexm1,4)];
HSVbinsNG = mexUpdateHist(HSVbinsNG,hsv',POI_NG,fadingRate,int32(wndSampling)); 

POI_G = [mseq.m1(groundIdx,3) mseq.m1(groundIdx,4)];
HSVbinsG = mexUpdateHist(HSVbinsG,hsv',POI_G,fadingRate,int32(wndSampling)); 

% G Points      
ntot = sum(HSVbinsG);                
log_Px = log((HSVbinsG+0.1)/(ntot+0.1*Nquan));     % Laplacian smoothing  

% NG points                        
ntot = sum(HSVbinsNG);                
log_PNGx = log((HSVbinsNG+1)/(ntot+Nquan));       % Laplacian smoothing             

llhngp = [];
llhp = [];
GMap = zeros([height width]);
for kj = (1+wndTesting):(height-wndTesting-1)
    for ki = (1+wndTesting):(width-wndTesting-1)
        % for every subarea
        llhng = 0;
        llh = 0;
        wi = kj;
        wj = ki;
    %   for wi = max([kj-wndTesting 1]):min([kj+wndTesting height])
    % for wj = max([ki-wndTesting 1]):min([ki+wndTesting width])
                bin_idx = hsv(wi,wj) + 1;       
                llhng = llhng + log_PNGx(bin_idx);   
                llh = llh + log_Px(bin_idx);                               
    %        end
    %   end

        llhp(kj,ki) = llh + log_p_g;
        llhngp(kj,ki) = llhng + log_p_ng;
        GMap(kj,ki) = llhp(kj,ki) - llhngp(kj,ki) ;

    end  
end

GndProb = exp(llhp)./(exp(llhp) + exp(llhngp));                                                  
%displayGndProb

end