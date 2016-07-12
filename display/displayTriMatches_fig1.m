% display matches 
figure(1),
for n=1:size(triMatch,1)     
    plot( [triMatch(n,3) triMatch(n,5)], [triMatch(n,4) triMatch(n,6)], '-','LineWidth',1,'Color', 'k');    
end
drawnow;
