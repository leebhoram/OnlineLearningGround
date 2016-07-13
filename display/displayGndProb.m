red = 2*max((GndProb-0.5),0);
green = 2*max((0.5-GndProb),0);

gmap = zeros([height width 3]);
gmap(2:end-1,2:end-1,1) = red;
gmap(2:end-1,2:end-1,2) = green;  

figure(4), hold off;
set(gcf,'Color','w'); 
imshow(gmap);

figure(5), hold off;
set(gcf,'Color','w'); 
imagesc(GMap);

