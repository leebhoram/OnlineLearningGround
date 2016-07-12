dbclear all;

mex('./mex/mexMatchFeatureB_Geiger.cpp', 'CXXFLAGS=-msse3 -fPIC -O');
mex('./mex/mexFindFeature_Geiger.cpp', '-O');
mex('./mex/mexPredictBackFlow.cpp', '-O');
mex('./mex/mexRGB2Gray218IndexedHSV.cpp', '-O');
mex('./mex/mexSolveArrow3.cpp', '-O');
mex('./mex/mexUpdateHist.cpp', '-O');
