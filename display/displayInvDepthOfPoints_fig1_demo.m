% display the depth in color
load depth_color.mat
depth = round((Zinv)*40);
depth(find(depth < 1)) = 1;
depth(find(depth > 64)) = 64;
rgb = c(depth,:);
figure(1),   
scatter( tri(zIdx,5), tri(zIdx,6), 12, rgb,'fill');
title('inverse depth (Red-close, Blue-far)');


