dbclear all;

mex('./mexMatchFeatureB_Geiger.cpp', 'CXXFLAGS=-msse3 -fPIC -O');
mex('./mexFindFeature_Geiger.cpp', '-O');
mex('./mexPredictBackFlow.cpp', '-O');
mex('./mexRGB2Gray218IndexedHSV.cpp', '-O');
mex('./mexSolveArrow3.cpp', '-O');
mex('./mexUpdateHist.cpp', '-O');
