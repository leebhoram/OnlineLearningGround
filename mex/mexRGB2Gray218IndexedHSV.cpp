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
#include "mex.h"
#include <stdint.h>
#include <math.h>

typedef unsigned char byte;

#define RGB2GRAY(r,g,b) (byte)round(0.299*(double)r + 0.587*(double)g + 0.114*(double)b)

//hsv[0]: 0<h<1
//hsv[1]: 0<s<1
//hsv[2]: 0<v<255
inline void RGB2HSV(byte r, byte g, byte b, float hsv[3]);

void convRGB2Gray218IndexedHSV(const byte *rgb,const int32_t* dims,byte* gray, byte* qHSV);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{      
   if (nrhs != 1 || nlhs != 2) 
        mexErrMsgTxt("Wrong number of arguments.");
    
    // read input arguments
    #define IM_IN           prhs[0]
    #define GRAY_OUT        plhs[0]
    #define INDHSV_OUT      plhs[1]
   
    unsigned char   *pIn;
    unsigned char   *pOut;
    unsigned char   *pOut2;
    int w, h;
   
    mxClassID classIn = mxGetClassID(IM_IN);
    
    if (classIn != mxUINT8_CLASS ) 
        mexErrMsgTxt("Invalid type for the input argument!");
   
    h = mxGetM(IM_IN);        
    w = mxGetN(IM_IN)/3;   
    pIn = (unsigned char *)mxGetData(IM_IN);  
    
           
    int32_t dims[3] = {h,w,3};
        // output argument
    GRAY_OUT =  mxCreateNumericArray(2,dims,mxUINT8_CLASS, mxREAL);
    pOut = (unsigned char *)mxGetData(GRAY_OUT);   
    
    INDHSV_OUT =  mxCreateNumericArray(2,dims,mxUINT8_CLASS, mxREAL);
    pOut2 = (unsigned char *)mxGetData(INDHSV_OUT);  
  
    convRGB2Gray218IndexedHSV(pIn,dims,pOut,pOut2);
    
    return;
}

void convRGB2Gray218IndexedHSV(const byte *rgb,const int32_t* dims,byte* gray, byte* qHSV)
{
    if (rgb == NULL || gray == NULL || qHSV == NULL)	return;

    int ch = dims[2]*dims[0];
    int h = dims[0];
    int w = dims[1];
    byte r,g,b;
    float hsv[3];
    byte q, q_;
    float vqlevel = 6.0/255.0;
	for (int32_t i =0; i<w;i++)
		for (int32_t j=0; j<h; j++){
			r = (rgb[j + i*h ]);
			g = (rgb[j + i*h + w*h]);
			b = (rgb[j + i*h + 2*w*h]);
			gray[i*h + j] = RGB2GRAY(r,g,b);

			RGB2HSV(r,g,b,hsv);
        
            q_ = ((byte)(vqlevel * hsv[2]));
            if (q_ == 6) q_--;
         
			q = (byte)(hsv[0]*6.0) + 6*((byte)(6.0*hsv[1])) + 36*q_;
            
			qHSV[i*h+j] = q;
		}

	return;
}


//hsv[0]: 0<h<1
//hsv[1]: 0<s<1
//hsv[2]: 0<v<255
inline void RGB2HSV(byte r, byte g, byte b, float hsv[3])
{
	byte maxval = 0,minval = 0, ind = 0;
	minval = g;
	maxval = r; ind = 1;
	if (r < g) {maxval = g; ind = 2; minval = r; };
	if (b > maxval) {maxval = b; ind = 3;}
	else if (b < minval) minval = b;

	hsv[2] = float(maxval);
	float diff = hsv[2] - (float)minval;

	// 0<s<1
	if (maxval == 0) hsv[1] = 0;
	else hsv[1] = (diff/hsv[2]);
    
    if (hsv[1] > 0.999) hsv[1] = 0.999;

	// 0<h<1
	switch(ind)
	{
	case 1:
		hsv[0] = (float)(g-b)/diff/6.0;
		if (hsv[0] < 0) {hsv[0] = hsv[0] + 5.0/6.0;}
		break;
	case 2: hsv[0] = (2.0/6.0 + (float)(b-r)/diff/6.0);
		break;
	case 3: hsv[0] = (4.0/6.0 + (float)(r-g)/diff/6.0);
		break;
	}

}

