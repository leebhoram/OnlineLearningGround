function tri = computeTriplet(mseq)

% getTripleMatch  
[Img2Idx, iprev, icur] = intersect( mseq.m1(:,6),mseq.m2(:,5));      
% Update age
deadIdx = setdiff(1:length(mseq.m1),iprev);
mseq.age(deadIdx) = 0;
mseq.age(iprev) = mseq.age(iprev) + 1;  
triMatch = double([mseq.m1(iprev,1:4)  mseq.m2(icur,3:4) iprev icur Img2Idx ...
            mseq.age(iprev) mseq.gndLLH(iprev) mseq.invDepthVar(iprev)]);   

% can be proved....
% check the 'straigint-ness' of matches
diff = ([triMatch(:,3:6)-triMatch(:,1:4) triMatch(:,5:6)-triMatch(:,1:2)]);
the1 = atan2(diff(:,1),diff(:,2));
the2 = atan2(diff(:,3),diff(:,4));
the12 = atan2(diff(:,5),diff(:,6));
a = [];
idxa = 0;
for k=1:size(diff,1)
    if norm(diff(k,5:6)) > 2                       
        if (abs(the1(k)-the2(k)) < 45*pi/180) &&  (abs(the2(k)-the12(k)) < 30*pi/180)  &&  (abs(the1(k)-the12(k)) < 30*pi/180)     
            idxa = idxa + 1;
            a(idxa,:) = [k -diff(k,6) diff(k,5) the12(k)]; %#ok<SAGROW>
        end
    else
        idxa = idxa + 1;
        a(idxa,:) = [k -diff(k,6) diff(k,5) the12(k)]; %#ok<SAGROW>
    end
end 
tri = triMatch(a(:,1),:);


