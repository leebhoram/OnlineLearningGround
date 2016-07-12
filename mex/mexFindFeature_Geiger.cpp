/* This is a modified matlab code from a part of libelas(c++).
 * The original code is published under GPL
 * and downloadable at http://www.cvlibs.net/software/libelas/
 * 
 * modified by Bhoram Lee (bhoram.lee@gmail.com)
 */
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <math.h>
#include "mex.h"

using namespace std;

const static int Scharr[3][3] = {{-3, 0, 3},{-10, 0, 10},{-3,0,3}};


const short int KB[5*5] = {-1, -1, -1, -1, -1,
											 -1,  1,  1,  1, -1,
											 -1,  1,  8,  1, -1,
											 -1,  1,  1,  1, -1,
											 -1, -1, -1, -1, -1};

const short int KC[5*5] = {-1, -1,  0,  1,  1,
											 -1, -1,  0,  1,  1,
											 0,  0,  0,  0,  0,
											 1,  1,  0, -1, -1,
											 1,  1,  0, -1, -1};

const int IdxFD[16][2] = { {5,  1},  {7,  1},
                     {1,  3},  {11, 3},
                     {3,  5},  {5,  5},
                     {7,  5},  {9,  5},
                     {3,  7},  {5,  7},
                     {7,  7},  {9,  7},
                     {1,  9},  {11,  9},
                     {5, 11},  {7,  11}};  // Feature descriptor index

typedef struct nmsPoints{
		short int u;   		// u-coordinate
		short int v;   		// v-coordinate
		short int val; 		// value
		short int c;        // class		
        unsigned char f[32];
		nmsPoints(){};
		nmsPoints(int u,int v,short int val,short int c):u(u),v(v),val(val),c(c){
           std::memset(f,0,32*sizeof(unsigned char));
        };        
}SnmsPoints;

int  findFeatures(const unsigned char* pin, int imw, int imh, std::vector<SnmsPoints>& pts, short int param[6]);
void scharr3x3_uint8Sat(const unsigned char* pin,int imw, int imh, char type, unsigned char* pout);
void convolve5x5Kernel(const unsigned char* pin, int imw, int imh, short int* pout, const short int K[5*5]);
void efficientNonMaxSupp_int16(const short int* pin, int imw, int imh, std::vector<SnmsPoints>& pts, short int param[4]);
void formFVec(const unsigned char* ph, const unsigned char* pv,  int imw, int imh, int x, int y, unsigned char fd[32]);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{ 
    //mexPrintf("int size:%d",sizeof(int));
   if (nrhs != 2 || nlhs != 2) 
        mexErrMsgTxt("Wrong number of arguments.");
    
    // read input arguments
    #define IMG_IN       prhs[0]
    #define PARAM        prhs[1]
    #define PTS_OUT      plhs[0]
    #define FTS_OUT      plhs[1]
   
    unsigned char *pIn;
    short int     *pOut1;
    unsigned char *pOut2;
    short int     *param;           
    char          type;           
    int w,h;

    h = mxGetN(IMG_IN);  // # of rows
    w = mxGetM(IMG_IN);  // # of cols
    
    //mexPrintf("M(h):%d, N(w):%d\n",h,w);
          
    mxClassID classIn = mxGetClassID(IMG_IN);
    mxClassID classType = mxGetClassID(PARAM);
    
    if (classIn != mxUINT8_CLASS) 
        mexErrMsgTxt("The first input argument must be a 2D uint8 array.");
    
    if (classType != mxINT16_CLASS) 
        mexErrMsgTxt("The second input argument must be a 1D INT16 array.");
    
    pIn = (unsigned char*)mxGetData(IMG_IN); // Get the pointer
    param = (short int *)mxGetData(PARAM);
    
    std::vector<SnmsPoints> vPts;
        
    int num = findFeatures(pIn, w, h, vPts, param);   
        
    if (num > 0)
    {     
        // output argument
        int dims[2] = {num,4};
        PTS_OUT =  mxCreateNumericArray(2,dims,mxINT16_CLASS, mxREAL);
        pOut1 = (short int *)mxGetData(PTS_OUT);   
              
        // output argument
        dims[0] = 32;
        dims[1] = num;
        FTS_OUT =  mxCreateNumericArray(2,dims,mxUINT8_CLASS, mxREAL);
        pOut2 = (unsigned char *)mxGetData(FTS_OUT);  
                        
        int idx = 0;
        for (vector<SnmsPoints>::iterator it=vPts.begin(); it!=vPts.end(); it++,idx++) {
            pOut1[idx+num*0 ] = (*it).u+1;
            pOut1[idx+num*1 ] = (*it).v+1;
            pOut1[idx+num*2 ] = (*it).val;
            pOut1[idx+num*3 ] = (*it).c;  

            std::memcpy(&pOut2[idx*32],(*it).f,32*sizeof(unsigned char));            
        }
        
    }
    else
    { 
        int dims[2] = {1,1};
        PTS_OUT =  mxCreateNumericArray(2,dims,mxINT16_CLASS, mxREAL);
        pOut1 = (short int *)mxGetData(PTS_OUT);    

        FTS_OUT =  mxCreateNumericArray(2,dims,mxUINT8_CLASS, mxREAL);
        pOut2 = (unsigned char *)mxGetData(FTS_OUT);  
        
        mexErrMsgTxt("Error while finding features");
    }
    
    vPts.clear();
    return;
}

int findFeatures(const unsigned char* pin, int imw, int imh, std::vector<SnmsPoints>& pts, short int param[6])
{
    if (pin == 0 || imw < 3 || imh < 3)
    return -1;
    
    unsigned char* gradv = new unsigned char [imw*imh];
    unsigned char* gradh = new unsigned char [imw*imh];
    scharr3x3_uint8Sat(pin,imw,imh,'v', gradv);
    scharr3x3_uint8Sat(pin,imw,imh,'h', gradh);
        
    if (param[0] > 0)
    {
        short int* B = new short int [imw*imh];
        convolve5x5Kernel(pin,imw,imh, B, KB);
        efficientNonMaxSupp_int16(B,imw,imh,pts,&(param[0]));
        delete [] B;
    }
    if (param[0+3] > 0)
    {
        short int* C = new short int [imw*imh];
        convolve5x5Kernel(pin,imw,imh, C, KC);
        efficientNonMaxSupp_int16(C,imw,imh,pts,&(param[0+3]));
        delete [] C;
    }
         
    int idx = 0;
    // Form feature vector
    if (pts.size() > 0)
    {       
        for (vector<SnmsPoints>::iterator it=pts.begin(); it!=pts.end(); it++, idx++) { 
            formFVec(gradh, gradv, imw, imh, (*it).u, (*it).v, (*it).f);
            //mexPrintf("%d) %d %d\n",idx+1, (*it).u, (*it).v);
        }
    }
  
    if (gradv) delete [] gradv;
    if (gradh) delete [] gradh;
  
    
    return idx;
}

void scharr3x3_uint8Sat(const unsigned char* pin, int imw, int imh, char type, unsigned char* pout)
{
	if (pin == 0 || pout == 0 || imw < 3 || imh < 3)
		return ;

	int exx = 1;
	int exy = 1;
	double scale = .1; // generally, should be double

	int exSize = (imw+2*exx)*(imh + 2* exy);

	unsigned char* exImg = new unsigned char [exSize];
	std::memset(exImg,0,exSize);

	for(int n=0; n<imh; n++)
	{
		memcpy((unsigned char*)(exImg + (n+exy)*(imw+2*exx)+exx),(unsigned char*)(pin+n*imw),imw);
	}

    if (type == 'v'){
        for (int j=0;j<imh;j++){
            for (int k=0;k<imw;k++){

                int conv = (Scharr[0][0] * exImg[(j+0)*(imw+2*exx) + k+0]
                           + Scharr[1][0] * exImg[(j+0)*(imw+2*exx) + k+1]
                           + Scharr[2][0] * exImg[(j+0)*(imw+2*exx) + k+2]
                           + Scharr[0][1] * exImg[(j+1)*(imw+2*exx) + k+0]
                           + Scharr[1][1] * exImg[(j+1)*(imw+2*exx) + k+1]
                           + Scharr[2][1] * exImg[(j+1)*(imw+2*exx) + k+2]
                           + Scharr[0][2] * exImg[(j+2)*(imw+2*exx) + k+0]
                           + Scharr[1][2] * exImg[(j+2)*(imw+2*exx) + k+1]
                           + Scharr[2][2] * exImg[(j+2)*(imw+2*exx) + k+2]);

                int val = round((double)conv * scale);

                if (val < -128) val = -128;
                else if (val > 127) val = 127;

                val += 128;
                pout[j*imw + k] = (unsigned char)val;
            }
        }
    }
    else
    {
        for (int j=0;j<imh;j++){
            for (int k=0;k<imw;k++){

                int conv = (Scharr[0][0] * exImg[(j+0)*(imw+2*exx) + k+0]
                           + Scharr[0][1] * exImg[(j+0)*(imw+2*exx) + k+1]
                           + Scharr[0][2] * exImg[(j+0)*(imw+2*exx) + k+2]
                           + Scharr[1][0] * exImg[(j+1)*(imw+2*exx) + k+0]
                           + Scharr[1][1] * exImg[(j+1)*(imw+2*exx) + k+1]
                           + Scharr[1][2] * exImg[(j+1)*(imw+2*exx) + k+2]
                           + Scharr[2][0] * exImg[(j+2)*(imw+2*exx) + k+0]
                           + Scharr[2][1] * exImg[(j+2)*(imw+2*exx) + k+1]
                           + Scharr[2][2] * exImg[(j+2)*(imw+2*exx) + k+2]);

                int val = round((double)conv * scale);

                if (val < -128) val = -128;
                else if (val > 127) val = 127;

                val += 128;
                pout[j*imw + k] = (unsigned char)val;
            }
        }
    }

    delete [] exImg;
	return;
}

void convolve5x5Kernel(const unsigned char* pin, int imw, int imh, short int* pout, const short int K[5*5])
{
	int exx = 2;
	int exy = 2;

	int exSize = (imw+2*exx)*(imh + 2* exy);

	unsigned char* exImg = new unsigned char [exSize];
	std::memset(exImg,0,exSize);

	for(int n=0; n<imh; n++)
	{
		memcpy((unsigned char*)(exImg + (n+exy)*(imw+2*exx)+exx),(unsigned char*)(pin+n*imw),imw);
	}

	int j=0;
	int k=0;
	for(j=0; j<imh; j++){
		for(k=0; k<imw; k++){
			short int val = 0;
			for(int y=0; y<5; y++){
				for(int x=0; x<5; x++){
					int idx = k+x;
					int idy = j+y;
					if (idx >= 0 && idy >= 0)
					{
						short int n1 = (short int) exImg[idy*(imw+2*2) + idx];
						val += n1 * K[y*5+x];
					}
				}
			}
			pout[j*imw + k] = val;
		}
	}

     delete [] exImg;
	return;
}


void efficientNonMaxSupp_int16(const short int* pin, int imw, int imh, std::vector<SnmsPoints>& pts, short int param[3])
{
	if (pin == 0 || param == 0)
		return ;

	short int nms_wnd = param[0];//3;
	short int tau = param[1]; //50;
	short int margin = 2*nms_wnd; //6;
	short int feature_class_flag_offset = param[2];
	int x = 0, y = 0;
    
	for (x = (nms_wnd + margin); x <= (imw - nms_wnd - margin); x +=(nms_wnd+1)){
		for (y = (nms_wnd + margin); y <= (imh - nms_wnd - margin); y +=(nms_wnd+1))	{

			int fminx = x;
			int fminy = y;
			int fmaxx = x;
			int fmaxy = y;
			short int fminval = pin[y*imw + x];
			short int fmaxval = fminval;

			for (int x1 = x; x1 <= (x + nms_wnd); x1++){
				for (int y1 = y; y1 <= (y + nms_wnd); y1++){
					short int curval = pin[y1*imw + x1];
                                      
					if(curval < fminval){
						fminx = x1;
						fminy = y1;
						fminval = curval;
					}
					else if (curval > fmaxval){
						fmaxx = x1;
						fmaxy = y1;
						fmaxval = curval;
					}
				}
			}

			bool failed = false;
			int x1 = 0;
			int y1 = 0;
			for (x1 = (fminx-nms_wnd); x1 <= std::min(fminx+nms_wnd,imw-margin-2); x1++){
				for (y1 = (fminy-nms_wnd); y1 <= std::min(fminy+nms_wnd,imh-margin-2); y1++){

					short int curval = pin[y1*imw + x1];
					if (curval < fminval) {
						if  ( (x1 < x) || (x1 > x + nms_wnd) || (y1 < y) || (y1 > y + nms_wnd)) {
							failed = true;
							break;
						}
					}
				}
				if (failed == true) break;
			}

			if (failed == false && fminval <= -tau)
			{
				//pout[y1*imw + x1] = (unsigned char)255;
				pts.push_back(SnmsPoints(fminx, fminy,fminval,0+feature_class_flag_offset));
			}

			failed = false;
			x1 = 0; y1 = 0;
			for (x1 = (fmaxx-nms_wnd); x1 <= std::min(fmaxx+nms_wnd,imw-margin-2); x1++){
				for (y1 = (fmaxy-nms_wnd); y1 <= std::min(fmaxy+nms_wnd,imh-margin-2); y1++){

					short int curval = pin[y1*imw + x1];

					if (curval > fmaxval) {
						if  ( (x1 < x) || (x1 > x + nms_wnd) || (y1 < y) || (y1 > y + nms_wnd)) {
							failed = true;
							break;
						}
					}
				}
				if (failed == true) break;
			}

			if (failed == false && fmaxval >= tau)
			{
				//pout[y1*imw + x1] = (unsigned char)255;
				pts.push_back(SnmsPoints(fmaxx,fmaxy,fmaxval,1+feature_class_flag_offset));
			}

			
		}
	}
    
	return;
}

void formFVec(const unsigned char* ph, const unsigned char* pv, int imw, int imh, int x, int y, unsigned char fd[32])
{
    
    for (int k=0;k<16;k++){
		int u = IdxFD[k][0] - 6;
		int v = IdxFD[k][1] - 6;

		fd[2*k] = ph[(y+v)*imw + (x+u)];
		fd[2*k+1] = pv[(y+v)*imw + (x+u)];
	}
	return;
}