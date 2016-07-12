function bstatic = is_static(mseq)

bstatic = false;

% test_static
abs_diff = abs(mseq.m2(:,3:4)-mseq.m2(:,1:2));
dx_lt1 = max(sum(abs_diff < 1))/length(mseq.m2);     
dx_lt2 = max(sum(abs_diff < 2))/length(mseq.m2);    
 
if dx_lt1 >= 0.8 && dx_lt2 > 0.9
   % disp('The scene is static.'); 
    bstatic = true;
else  % ==== when V > 0   
    bstatic = false;
end         
