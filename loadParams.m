% [Global Group-1] Camera or image properties
global width
global height
global f cv cu h
global dt;

dt = 0.1;  % 10Hz 
im = imread(imagelist{1});  
width = 621;
height = 188;   
h = 1.68;
f = 718.856; cu = 607.1928; cv = 185.2157; 
f = f/2; cu = cu/2; cv = cv/2;      

% to keep track of matches 
MatchSequence = struct('m1',[],'m2',[],'gndLLH',[],'invDepth',[],...
                       'invDepthPrev',[],'invDepthFlag',[],'invDepthPred',[],'invDepthEst',[],'invDepthVar',[]);

% RANSAC Param 
NUM = 3; % check NUM_ITER
NUM_ITER = ceil(log(1-0.95)/log(1-(1-0.65)^NUM));  
global Ransac_param;
Ransac_param= struct('req_numIter', NUM_ITER,...
                      'req_numData',0.35,...                     % 'sig_dist',0.015,...
                      'sig_dist',0.005);
                     % 'req_minErr',0.05);
rng('shuffle');  
clear NUM NUM_ITER

% statistics of color
global Nquan % level of quantization
CHquan = 6;
Nquan = CHquan^3;
HSVbinsNG = zeros(1,Nquan);
HSVbinsG = zeros(1,Nquan);
GndProb = [];

% groun truth data
loadRef % for comparision, also used to set initial motion


% filter for velocity
Vhat = norm(ref_path(startIdx+1,:)-ref_path(startIdx,:))*10;
VF = KF1D;
VF = VF.initialize(Vhat,0.15,3);

% filter for angular velocity
alphaHat = -(ref_eulr(startIdx+1,3) - ref_eulr(startIdx,3))*10;
wF = KF1D;
wF = wF.initialize(alphaHat,2,5);
alpha_prev = alphaHat;

Phat_mat = zeros(1,3);
yawHat = 0;
P = [0 0];



