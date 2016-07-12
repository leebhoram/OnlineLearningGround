if fn > (startIdx + 1)
    [y,id] = sortrows(triMatch(:,7:8),2);
    temp2 = zeros(length(MatchSequence.m1),1);

    temp2(y(:,2)) = MatchSequence.age(y(:,1));
    MatchSequence.age = temp2;

    temp2(y(:,2)) = MatchSequence.gndLLH(y(:,1));
    MatchSequence.gndLLH = temp2;
    
%     temp3 = 0.5*ones(length(MatchSequence.m1),1);
%     temp3(y(:,2)) = MatchSequence.gndDist(y(:,1));
%     MatchSequence.gndDist = temp3;
%         
%     temp2(y(:,2)) = MatchSequence.gndLLH2(y(:,1));
%     MatchSequence.gndLLH2 = temp2;

%     temp2(y(:,2),:) = MatchSequence.prevStrongGnd(y(:,1),:);
%     MatchSequence.prevStrongGnd = temp2;

    temp2(y(:,2)) = MatchSequence.invDepthFlag(y(:,1));
    MatchSequence.invDepthFlag = temp2;           

    temp2(y(:,2)) = MatchSequence.invDepth(y(:,1));
    MatchSequence.invDepth = temp2;   
    MatchSequence.invDepthPrev = temp2; 

    temp2(y(:,2)) = MatchSequence.invDepthEst(y(:,1));
    MatchSequence.invDepthEst = temp2;
    
    temp2(y(:,2)) = MatchSequence.invDepthPred(y(:,1));
    MatchSequence.invDepthPred = temp2;

    temp3 = 2*ones(length(MatchSequence.m1),1);
    temp3(y(:,2)) = MatchSequence.invDepthVar(y(:,1));
    MatchSequence.invDepthVar = temp3;

end
