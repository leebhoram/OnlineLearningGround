/* This is a modified matlab code from a part of libviso2(c++).
 * The original code is published under GPL
 * and downloadable at http://www.cvlibs.net/software/libviso/
 * 
 * modified by Bhoram Lee (bhoram.lee@gmail.com)
 */
#include <stdlib.h>
#include <stdint.h>
#include <cstring>
#include <vector>
#include <math.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "mex.h"

using namespace std;
                     
typedef struct Match{
	short int idx1;
	short int x1;
	short int y1;
	short int idx2;
	short int x2;
	short int y2;
        
    short int sWnd[4];
    short int cost;
    
	Match(short int i1, short int x1, short int y1, short int i2, short int x2, short int y2):idx1(i1), x1(x1), y1(y1), idx2(i2), x2(x2),y2(y2) 
    {
        sWnd[0] = 0;
        sWnd[1] = 1;
        sWnd[2] = 2;
        sWnd[3] = 3;
    };
}SMatch;
                     
typedef struct nmsPoints{
		short int u;   		// u-coordinate
		short int v;   		// v-coordinate
		short int val; 		// value
		short int c;        // class     
        unsigned char* f;
        
        // matching info
		short int prev_idx;
		short int prev_u;
		short int prev_v;
		short int next_idx;
		short int next_u;
		short int next_v;
		int score;        
        
		nmsPoints(){};
		nmsPoints(int u,int v,short int val,short int c):u(u),v(v),val(val),c(c){
           prev_idx = -1; next_idx = -1;score = 0;
          // f = (short int*)_mm_malloc(32*sizeof(unsigned char),16);
           posix_memalign((void**)&f, 32,  32 * sizeof(unsigned char));
        };        
}SnmsPoints;

typedef union
{
    __m128i v;
    unsigned char a8[16];
    unsigned short int a16[8];
    uint32_t a32[4];
} U128;

int matchFeatureB(std::vector<SnmsPoints>& pts1,std::vector<SnmsPoints>& pts2, std::vector<SMatch> &match, short int* searchRect, double param[13]);
void createBinIndex(std::vector<SnmsPoints>& p,vector<int> *k,double param[3]);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{ 
    //mexPrintf("in:%d out:%d",nrhs, nlhs);
    
   if (nrhs != 6 || nlhs != 1) 
        mexErrMsgTxt("Wrong number of arguments.");
    
    // read input arguments
    #define P1_IN       prhs[0]
    #define P2_IN       prhs[1]
    #define F1_IN       prhs[2]
    #define F2_IN       prhs[3]
    #define SREGION        prhs[4]
    #define PARAM         prhs[5]
    #define M_OUT         plhs[0]
   
    short int       *p1In;
    short int       *p2In;
    unsigned char   *f1In;
    unsigned char   *f2In;
    short int       *srectIn;
    short int       *pOut1;
    short int       *pOut2;
    double          *param;      
    int              type;  
    int n1, n2;
   
    mxClassID classP1In = mxGetClassID(P1_IN);
    mxClassID classP2In = mxGetClassID(P2_IN);
    mxClassID classSRIn = mxGetClassID(SREGION);
    
    mxClassID classF1In = mxGetClassID(F1_IN);
    mxClassID classF2In = mxGetClassID(F2_IN);
    
    mxClassID classParam = mxGetClassID(PARAM);
    
    if (classSRIn != mxINT16_CLASS || classP1In != mxINT16_CLASS || classP2In != mxINT16_CLASS || classF1In != mxUINT8_CLASS || classF2In != mxUINT8_CLASS ) 
        mexErrMsgTxt("Invalid type for the input arguments!");
    
    if (classParam != mxDOUBLE_CLASS) 
        mexErrMsgTxt("The second input argument must be a 1D INT16 array.");
    
    n1 = mxGetM(P1_IN);       // # of points    
    if (mxGetN(P1_IN) != 4)   mexErrMsgTxt("A point must include two elements.");// 
    
    n2 = mxGetM(P2_IN);       // # of points
    if (mxGetN(P2_IN) != 4)   mexErrMsgTxt("A point must include two elements.");//   
   
    //mexPrintf("%d %d\n",n1, mxGetM(SREGION));
    
    if (mxGetM(SREGION) != n1)   mexErrMsgTxt("The region must include the same number of element to the first point set.");// 
    if (mxGetN(SREGION) != 4)    mexErrMsgTxt("The region must include four elements per row.");// 
    
    p1In = (short int*)mxGetData(P1_IN);  
    p2In = (short int*)mxGetData(P2_IN);  
    
    f1In = (unsigned char*)mxGetData(F1_IN);  
    f2In = (unsigned char*)mxGetData(F2_IN);  
    
    srectIn = (short int*)mxGetData(SREGION);  
    
    
    if (mxGetN(F1_IN) != n1 || mxGetN(F2_IN) != n2)   
        mexErrMsgTxt("The length of the points array and the features array must be the same.");// 
    
    param = (double *)mxGetData(PARAM);
    
    std::vector<SnmsPoints> vPts1;
    std::vector<SnmsPoints> vPts2;    
    std::vector<SMatch>     vM; 
    
   
    int numMatch = 0;
   
    int n=0 ;
    for (n=0; n<n1; n++) {
        vPts1.push_back(SnmsPoints(p1In[n+n1*0]-1,p1In[n+n1*1]-1,p1In[n+n1*2],p1In[n+n1*3]));        
    }
    //mexPrintf("n1:%d\n",vPts1.size());
    
    n=0;
    for ( vector<SnmsPoints>::iterator it=vPts1.begin(); it!=vPts1.end(); it++, n++) {
        for (int idx = 0; idx < 32; idx++){
            (*it).f[idx] = f1In[n*32 + idx];            
        }        
    }
    
    for (n=0; n<n2; n++) {
        vPts2.push_back(SnmsPoints(p2In[n+n2*0]-1,p2In[n+n2*1]-1,p2In[n+n2*2],p2In[n+n2*3]));
    }  
    //mexPrintf("n2:%d\n",vPts2.size());
    
    n=0;   
    for ( vector<SnmsPoints>::iterator it2=vPts2.begin(); it2!=vPts2.end(); it2++, n++) {
          for (int idx = 0; idx < 32; idx++){
             (*it2).f[idx] = f2In[n*32 + idx]; 
          } 
    }

    if (n1 > 0 && n2 > 0)  numMatch = matchFeatureB(vPts1, vPts2, vM, srectIn, param); 
    else                   numMatch = 0;
        
    if (numMatch > 0)
    {
        // mexPrintf("Number of Matches: %d\n",vM.size());
    
        // output argument
        int dims[2] = {numMatch,6};
        M_OUT =  mxCreateNumericArray(2,dims,mxINT16_CLASS, mxREAL);
        pOut1 = (short int *)mxGetData(M_OUT);   

        
        int idx = 0;
        for (vector<SMatch>::iterator it=vM.begin(); it!=vM.end(); it++,idx++) {
                      
           pOut1[idx+numMatch*0 ] = (*it).idx1+1;         
           pOut1[idx+numMatch*1 ] = (*it).idx2+1;
           pOut1[idx+numMatch*2 ] = (*it).x1+1;
           pOut1[idx+numMatch*3 ] = (*it).y1+1;
           pOut1[idx+numMatch*4 ] = (*it).x2+1;
           pOut1[idx+numMatch*5 ] = (*it).y2+1;       
        }
         
    }        
    else
    { 
        int dims[2] = {1,1};
        M_OUT =  mxCreateNumericArray(2,dims,mxINT16_CLASS, mxREAL);
        pOut1 = (short int *)mxGetData(M_OUT);    
        
        mexErrMsgTxt("Error while matching features");
    }
    
    for ( vector<SnmsPoints>::iterator it=vPts1.begin(); it!=vPts1.end(); it++) {
        free((*it).f);   
    }
    
    for ( vector<SnmsPoints>::iterator it=vPts2.begin(); it!=vPts2.end(); it++) {
        free((*it).f);   
    }
       
    vPts1.clear();
    vPts2.clear();
    vM.clear();
    
    return;
}

int matchFeatureB(std::vector<SnmsPoints>& pts1,std::vector<SnmsPoints>& pts2, std::vector<SMatch> &match, short int* searchRect, double param[13])
{
    double Vz = param[0];
    double w = param[1];
    int width = (int)param[2];
    int height = (int)param[3];
    /*
    double dt = param[4];
    double focl = param[5];  
     **/       
    double cu = param[6];
    double cv = param[7];
    
    
    int binSize =   (int)param[10];
    int nxbin =   (int)param[11];
    int nybin =   (int)param[12];
    int thre = (int)param[9];
    
	int wnd_minx = 1;
	int wnd_maxx = -1;
	int wnd_miny = 1;
	int wnd_maxy = -1;
    int R = (int)param[8];
     
    __m128i xmm1,xmm2;
    U128 xmm3,xmm4;
    
	int nBin = (int32_t)(4*nxbin*nybin);                  
	vector<int> *vp2 = new vector<int>[nBin];
	createBinIndex(pts2,vp2,&param[10]);
    
           
    // forward matching only
    int N = pts1.size();
	int n1 = 0;
	for (vector<SnmsPoints>::iterator it=pts1.begin(); it!=pts1.end(); it++, n1++) {
		int x = (*it).u;
		int y = (*it).v;
		int c = (*it).c;

		float u = (float)x - (float)cu;
		float v = (float)y - (float)cv;
    
		xmm1 = _mm_load_si128((__m128i*)((*it).f));
		xmm2 = _mm_load_si128((__m128i*)((*it).f+16));
        
		int us = 0, uf = 0, vs = 0, vf = 0;
        
        wnd_minx = searchRect[N*0 + n1];
        wnd_maxx = searchRect[N*1 + n1];
        wnd_miny = searchRect[N*2 + n1];
        wnd_maxy = searchRect[N*3 + n1]; 
            
        if ( (wnd_minx + wnd_maxx + wnd_miny + wnd_maxy) > 0)
        {
            us = std::max(wnd_minx,0);
            uf = std::min(wnd_maxx,width-1);
            vs = std::max(wnd_miny,0);
            vf = std::min(wnd_maxy,height-1);

            // bins of interest
            int u_bin_min = min(max((int)floor(us/(float)binSize),0),nxbin-1);
            int u_bin_max = min(max((int)floor(uf/(float)binSize),0),nxbin-1);
            int v_bin_min = min(max((int)floor(vs/(float)binSize),0),nybin-1);
            int v_bin_max = min(max((int)floor(vf/(float)binSize),0),nybin-1);

            int minval = 100000;
            int next_idx = -1;
            for (int ubinIdx = u_bin_min; ubinIdx <= u_bin_max; ubinIdx++){
                for (int vbinIdx = v_bin_min; vbinIdx <= v_bin_max; vbinIdx++){
                    int p2_bin_ind = (c*nybin+vbinIdx)*nxbin+ubinIdx;

                    for (vector<int>::iterator it2=vp2[p2_bin_ind].begin(); it2!=vp2[p2_bin_ind].end(); it2++) {
                        int x2 = pts2[*it2].u;
                        int y2 = pts2[*it2].v;

                        if (x2 >= us && x2 <= uf && y2 >= vs && y2 <= vf)
                        {
                           int32_t sad = 0;

                            xmm3.v = _mm_load_si128((__m128i*)(pts2[*it2].f));
                            xmm4.v = _mm_load_si128((__m128i*)(pts2[*it2].f+16));
                            xmm3.v = _mm_sad_epu8(xmm1,xmm3.v);
                            xmm4.v = _mm_sad_epu8(xmm2,xmm4.v);
                            sad = xmm3.a16[0]+xmm3.a16[4]+xmm4.a16[0]+xmm4.a16[4];
                            

                            if (sad < minval && sad < thre){
                                minval = (int)sad;
                                next_idx = (*it2) ;
                            }
                        } // end of if: Range Check
                    } // end of for
                } // end of for : v_bin
            } // end of for : u
            
          if (next_idx >= 0 && minval < thre){
                (*it).next_idx = next_idx;
                (*it).next_u = pts2[next_idx].u;
                (*it).next_v = pts2[next_idx].v;

                if (pts2[next_idx].prev_idx < 0) {
                    pts2[next_idx].prev_idx = (short int)n1;
                    pts2[next_idx].prev_u = (short int)x;
                    pts2[next_idx].prev_v = (short int)y;
                    pts2[next_idx].score = minval;
                }
                else
                {
                    if (minval < pts2[next_idx].score)
                    {
                        pts2[next_idx].prev_idx = (short int)n1;
                        pts2[next_idx].prev_u = (short int)x;
                        pts2[next_idx].prev_v = (short int)y;
                        pts2[next_idx].score = minval;
                    }
                }
            }
        }
	} // end of for: it1

      
	int n2 = 0;
    int num = 0;
	for (vector<SnmsPoints>::iterator it=pts2.begin(); it!=pts2.end(); it++, n2++) {
		if ((*it).prev_idx > 0) {
			match.push_back(SMatch((*it).prev_idx,(*it).prev_u,(*it).prev_v,n2,(*it).u,(*it).v));
            num++;
        }
	}

    for (int n=0;n<nBin;n++) {
		vp2[n].clear();
    }
    
	delete [] vp2;
    
    return num;
}

// param[0] = binSize
// param[1] = number of x bins
// param[2] = number of y bins
void createBinIndex(std::vector<SnmsPoints>& p,vector<int> *k, double param[3])
{
    int n1 = 0;
    int binSize = (int)param[0];
    int nxbin =   (int)param[1];
    int nybin =   (int)param[2];
  
    
	for (vector<SnmsPoints>::iterator it=p.begin(); it!=p.end(); it++, n1++) {
		int u = (*it).u;
		int v = (*it).v;
		int c = (*it).c;
        
		// compute row and column of bin to which this observation belongs
		int u_bin = min((int)floor((float)u/(float)binSize),nxbin-1);
		int v_bin = min((int)floor((float)v/(float)binSize),nybin-1);

		// save index
		k[(c*nybin+v_bin)*nxbin+u_bin].push_back(n1);   
        
    }    
      
    return;
}
