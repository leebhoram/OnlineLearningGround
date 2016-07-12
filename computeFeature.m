function [gray,hsv, pts,fvec] = computeFeature(img)

% Detection Parameters as in libviso
nms_wb = 3; % negative if not in use
nms_wc = 3;
nms_tau =50;
FDParam = int16([nms_wb nms_tau 2 nms_wc nms_tau 0]);

[gray, hsv] = mexRGB2Gray218IndexedHSV(img);    
[pts, fvec] =  mexFindFeature_Geiger(gray',FDParam); % pts: [1] x [2] y [3] nms val [4] class (0 or 1)

end