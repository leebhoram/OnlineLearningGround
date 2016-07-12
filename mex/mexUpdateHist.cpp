// Copyright (c) 2015 Bhoram Lee
// 
// This file is part of OLG.
// Authors: Bhoram Lee
// 
// OLG is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
// updatedHist = mexUpdateHist(Hist,data,POI,decayingRate,nWnd); 
// 
// INPUT
// Hist
//   1D double array
// nBin
//   number of bins of the histogram (positive integer (uint8))
// data
//   2D (quantized image) data for builing to the histogram
// POI
//   points of interest (2 x [# of points])
// decayingRate
//   the decaying rate for old histogram (0-1, double)
// nWnd
//   the window size (half wnd) for collecting data (positive integer)
//
// OUTPUT
// updatedHist 
//   1D double array

#include <string.h>
#include "mex.h"

using namespace std;
                     
void updateHist(const double* pHist,double* pHistNew,int nBin,unsigned char* data, int w, int h, short int *arrPOIx, short int *arrPOIy, int nPOI, double decayingRate, int nNeighborWnd);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{ 
    
   if (nrhs != 5 || nlhs != 1) 
        mexErrMsgTxt("Wrong number of arguments.");
    
    // read input arguments
    #define HIST_IN       prhs[0]
    #define DATA_IN       prhs[1]
    #define PTS_IN        prhs[2]
    #define DECAYRATE_IN  prhs[3]
    #define WND_IN        prhs[4]
    #define HIST_OUT      plhs[0]
   
    double          *phistIn;
    unsigned char   *pdataIn;
    short int       *pPtsIn;
    double          *phistOut;   
    int             dim1, dim2, nBin, nWnd=0;
    double          decayRate=0.0;
   
    mxClassID classHistIn = mxGetClassID(HIST_IN);
    mxClassID classDataIn = mxGetClassID(DATA_IN);
    mxClassID classPtsIn = mxGetClassID(PTS_IN);
        
    if (classHistIn != mxDOUBLE_CLASS || classDataIn != mxUINT8_CLASS || classPtsIn != mxINT16_CLASS ) 
        mexErrMsgTxt("Invalid type for the input arguments!");
          
    // Check dimensions
    nBin = mxGetN(HIST_IN);
    if (nBin < 2)
        mexErrMsgTxt("The first input must be a 1xNUM_BIN array.");
    
    dim1 = mxGetM(DATA_IN);       // row 
    dim2 = mxGetN(DATA_IN);       // dim
    
    int nPOI = 0;
    if (mxGetN(PTS_IN) == 2)   
    {
        nPOI = mxGetM(PTS_IN);
    }
    else
    {
        mexErrMsgTxt("The third input must be a Nx2 matrix.");// 
    }
   
    // Input Arguments
    phistIn = mxGetPr(HIST_IN);  
    pdataIn = (unsigned char*)mxGetData(DATA_IN);  
    pPtsIn = (short int*)mxGetData(PTS_IN);  
   
    nWnd = *(int*)mxGetData(WND_IN);  
    decayRate = *(double*)mxGetData(DECAYRATE_IN);      
        
    // Get the Output Pointer
    HIST_OUT =  mxCreateDoubleMatrix(1,nBin,mxREAL);
    phistOut = (double *)mxGetPr(HIST_OUT);   
        
    updateHist(phistIn, phistOut, nBin, pdataIn, dim1, dim2, pPtsIn, (short int*)&(pPtsIn[nPOI]), nPOI, decayRate, nWnd);
     
    return;
}

void updateHist(const double* pHist,double* pHistNew, int nBin,unsigned char* data, int w, int h, short int *arrPOIx, short int *arrPOIy, int nPOI, double decayingRate, int nNeighborWnd)
{
	if (!pHist || !pHistNew || !data || !arrPOIx || !arrPOIy || nBin < 2 || nPOI < 0 || nNeighborWnd < 0)
		return;

	double N = (double)(nBin*(2*nNeighborWnd+1)*(2*nNeighborWnd+1));
	double *tempHist = new double [nBin];
    memset(tempHist,0,nBin*sizeof(double));

	for (int i = 0; i < nPOI; i++)
	{
		for (int j=(arrPOIy[i] - nNeighborWnd-1);j<=(arrPOIy[i] + nNeighborWnd-1);j++)
			for (int k=(arrPOIx[i] - nNeighborWnd-1);k<=(arrPOIx[i] + nNeighborWnd-1);k++){
				if (j >=0 && j < h && k >=0 && k < w){
					int b = (int)data[j*w + k];                   
				
					if ( b>=0 && b< nBin)
					{
						tempHist[b] = tempHist[b]+1.0;
					}
				}
			}
	}

	for (int b=0;b<nBin;b++)
	{
		pHistNew[b] = (decayingRate*pHist[b] + tempHist[b]);
	}
	
	delete [] tempHist;

	return;
}



