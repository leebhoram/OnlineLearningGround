% display the slope (scale estimation)
figure(3);  hold off;
set(gcf,'Color','w');
plot(lowerHalf(:,1),lowerHalf(:,2),'.','Color',[0.5 0.5 0.5]); hold on;
scatter( lowerHalf(inliers,1), lowerHalf(inliers,2),'bo'); 
xlabel('normalized y coordinate');
ylabel('inverse depth meas (multiplied by V)');


