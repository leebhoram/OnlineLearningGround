function mseq = updateInvDepth(mseq)
       
ZinvR = 0.6;   
priorInvDepthIdx = find((mseq.invDepthVar < 5) & (mseq.invDepthFlag == 1)) ;   
noPriorInvDepthIdx = setdiff(1:length(mseq.invDepthFlag),priorInvDepthIdx);

if ~isempty(noPriorInvDepthIdx)
    mseq.invDepthEst(noPriorInvDepthIdx) = mseq.invDepth(noPriorInvDepthIdx);
    mseq.invDepthVar(noPriorInvDepthIdx) = 0.95*ones(length(noPriorInvDepthIdx),1);
end

if ~isempty(priorInvDepthIdx) 
    R_plus_P = ZinvR + mseq.invDepthVar(priorInvDepthIdx);
    mseq.invDepthEst(priorInvDepthIdx) = 1./R_plus_P.*(ZinvR*mseq.invDepthPred(priorInvDepthIdx) + ...
                                                  mseq.invDepthVar(priorInvDepthIdx).*mseq.invDepth(priorInvDepthIdx));
    mseq.invDepthVar(priorInvDepthIdx) = ZinvR*mseq.invDepthVar(priorInvDepthIdx)./R_plus_P;
end  