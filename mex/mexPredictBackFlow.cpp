#include <string.h>
#include <cmath> 
#include "matrix.h"
#include "mex.h"

using namespace std;
                    
void predict(const double* px0, int nx, const double* invDepth0, double* xRect,double* xdot, double* types, double param[9]);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{     
   if (nrhs != 3 || nlhs != 3) 
        mexErrMsgTxt("Wrong number of arguments.");
   
   // read input arguments
    #define X0_IN        prhs[0]
    #define INVD_IN      prhs[1]
    #define PARAM_IN     prhs[2]
    #define XRECT_OUT    plhs[0]
    #define XDOT_OUT     plhs[1]
    #define TYPE_OUT     plhs[2]
   
    double          *px0;
    double          *pInvD;
    int             dim1, dim2;
    double          *pParamIn;
    double          *pTypeOut;
    double          *pXRectOut;
    double          *pXdotOut;
   
    mxClassID cX0 = mxGetClassID(X0_IN);
    mxClassID cInvD = mxGetClassID(INVD_IN);
    mxClassID cParamIn = mxGetClassID(PARAM_IN);
        
    if (cX0 != mxDOUBLE_CLASS || cInvD != mxDOUBLE_CLASS ) 
        mexErrMsgTxt("Invalid type for the input arguments!");
                
    dim1 = mxGetN(X0_IN);       // row 
    dim2 = mxGetM(X0_IN);       // dim
    
              
   // Input Arguments
    if (dim1 != mxGetN(INVD_IN))
        pInvD = 0;
    else    
        pInvD = mxGetPr(INVD_IN);
    
    px0 = mxGetPr(X0_IN);    
    pParamIn = mxGetPr(PARAM_IN); 
     
    // Get the Output Pointer
    XRECT_OUT =  mxCreateDoubleMatrix(dim1,4,mxREAL);
    pXRectOut = (double *)mxGetPr(XRECT_OUT);  
    XDOT_OUT =  mxCreateDoubleMatrix(dim1,2,mxREAL);
    pXdotOut = (double *)mxGetPr(XDOT_OUT); 
    TYPE_OUT =  mxCreateDoubleMatrix(dim1,1,mxREAL);
    pTypeOut = (double *)mxGetPr(TYPE_OUT); 
   
    predict(px0, dim1, pInvD, pXRectOut, pXdotOut, pTypeOut, pParamIn);

    return;
}

void predict(const double* px0, int nx, const double* invDepth0, double* xRect,double* xdot, double* types, double param[9])
{
    if (px0 == 0 || nx < 1 ||  xRect == 0 || xdot == 0 || types == 0)
        return;
    
    double f =  param[0];  // focal length
    double cv = param[1];
    double cu = param[2]; 
    double dt = param[3];
    double width = param[4];
    double height = param[5];
    double R = param[6];
    double w = param[7];
    double V = param[8];
    
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    
    if (fabs(w) > 0.1)  R = 2*R;
    
    for (int k=0;k<nx;k++){   
        double u = (px0[2*k+0] - cu)/f;
        double v = (px0[2*k+1] - cv)/f;
        
        if (invDepth0 > 0 && invDepth0[k] > 0 && V > 0){

            xdot[k+0] =  u*invDepth0[k]*V - (1 + u*u)*w;
            xdot[k+1*nx] =  v*invDepth0[k]*V - (u*v)*w;

            double x1 = px0[2*k+0] - f*xdot[k]*dt;
            double x2 = px0[2*k+1] - f*xdot[k+1*nx]*dt;

            if (x1 < 1 || x2 < 1 || x1 > width || x2 > height){
                xRect[k] = 0;
                xRect[1*nx + k] = 0;
                xRect[2*nx + k] = 0;
                xRect[3*nx + k] = 0;
                types[k]   = -1;
            }
            else {
                xmin = fmax((x1-2*R), 1.0f);
                xmax = fmin((x1+2*R), width);
                ymin = fmax((x2-2*R), 1.0f);
                ymax = fmin((x2+2*R), height);

                xRect[k] = xmin;
                xRect[1*nx + k] = xmax;
                xRect[2*nx + k] = ymin;
                xRect[3*nx + k] = ymax;
                types[k] = 0;
            }
        }
        else
        {
            double xdot_v01 =  0 - (1 + u*u)*w;
            double xdot_v02 =  0 - (u*v)*w;

            double xdot_v11 =  V/5*u + xdot_v01;
            double xdot_v12 =  V/5*v + xdot_v02;

            xdot[k] = -f*xdot_v11*dt;
            xdot[k+1*nx] = -f*xdot_v12*dt;

            double x11 = px0[2*k+0] - f*xdot_v01*dt;
            double x12 = px0[2*k+1] - f*xdot_v02*dt;

            double x21 = px0[2*k+0] + xdot[k];
            double x22 = px0[2*k+1] + xdot[k+1*nx];

            xmin = fmin(x11, x21);
            xmax = fmax(x11, x21);
            ymin = fmin(x12, x22);
            ymax = fmax(x12, x22);

            if  (xmin < 1 || ymin < 1|| xmax > width || ymax > height){
                xRect[0*nx + k] = 0;
                xRect[1*nx + k] = 0;
                xRect[2*nx + k] = 0;
                xRect[3*nx + k] = 0;
                types[k]   = -1;
            }
            else{      

                xmin = fmax((xmin-R), 1.0f);
                xmax = fmin((xmax+R), width);
                ymin = fmax((ymin-R), 1.0f);
                ymax = fmin((ymax+R), height);

                xRect[0*nx + k] = xmin;
                xRect[1*nx + k] = xmax;
                xRect[2*nx + k] = ymin;
                xRect[3*nx + k] = ymax;
                types[k] =1;
            }
        }        
    }
 
    return;
}

