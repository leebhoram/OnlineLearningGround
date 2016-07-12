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
// mexSolveArrow 
// 
// Solve Ax = b where A has the upper-left directed arraow structure
// 
//  * * * * *
//  * * 0 0 0
//  * 0 * 0 0
//  * 0 0 * 0
//  * 0 0 0 *

#include <string.h>
#include "matrix.h"
#include "mex.h"

using namespace std;
                    
void solveArrowLU3(const double* pA11, const double* pA12, const double* pA22, int n1, int n2, const double* b, double* x);
void inverse3x3Matrix(double *A_i, const double *A_);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{     
   if (nrhs != 4 || nlhs != 1) 
        mexErrMsgTxt("Wrong number of arguments.");
    
    // read input arguments
    #define A11_IN       prhs[0]
    #define A12_IN       prhs[1]
    #define A22_IN       prhs[2]
    #define B_IN       prhs[3]
    #define X_OUT      plhs[0]
   
    double          *pA11In;
    double          *pA12In;
    double          *pA22In;
    double          *pbIn;
    int             dim1, dim2;
    double          *pxOut;
   
    mxClassID classA11In = mxGetClassID(A11_IN);
    mxClassID classA12In = mxGetClassID(A12_IN);
    mxClassID classA22In = mxGetClassID(A22_IN);
    mxClassID classBIn = mxGetClassID(B_IN);
        
    if (classA11In != mxDOUBLE_CLASS || classA12In != mxDOUBLE_CLASS || classA22In != mxDOUBLE_CLASS || classBIn != mxDOUBLE_CLASS ) 
        mexErrMsgTxt("Invalid type for the input arguments!");
                
    dim1 = mxGetN(A11_IN);       // row 
    dim2 = mxGetM(A22_IN);       // dim
          
    // Input Arguments
    pA11In = mxGetPr(A11_IN);  
    pA12In = mxGetPr(A12_IN);  
    pA22In = mxGetPr(A22_IN);  
    pbIn = mxGetPr(B_IN); 
         
    // Get the Output Pointer
    X_OUT =  mxCreateDoubleMatrix(dim1+dim2,1,mxREAL);
    pxOut = (double *)mxGetPr(X_OUT);   
        
    solveArrowLU3(pA11In,pA12In,pA22In,dim1,dim2, pbIn, pxOut);
     
    return;
}


void solveArrowLU3(const double* pA11, const double* pA12, const double* pA22, int n1, int n2, const double* b, double *x)
{
	if (pA11 == 0 || pA12 == 0 || pA22 == 0 || b == 0 || n1 != 3 || n2 < 5)
		return;

	// A11 = 3 x 3 
    // A12 = 3 x [n2]
	
	double A12u[3] = {0, 0, 0};
	
	double A12v[3*3] = {0,0,0, 0,0,0, 0,0,0};
    double A_[3*3] =  {0,0,0, 0,0,0, 0,0,0} ; // A11 - A12v 
     double A_inv[3*3] =  {0,0,0, 0,0,0, 0,0,0} ;   // inv(A11 - A12v) 
     double b_[3] =  {0, 0, 0};
     double *u = new double[n2];
    double *v = new double[n2*3];
	
    for (int i=0;i<n2;i++)
	{
      
	    u[i] = b[i+n1]/pA22[i];
	    v[i*3+0] = pA12[i*3+0]/pA22[i];
        v[i*3+1] = pA12[i*3+1]/pA22[i];
        v[i*3+2] = pA12[i*3+2]/pA22[i];
        
       
        // A12*u  [3x1]
	    A12u[0] = A12u[0] + pA12[i*3+0]*u[i];
        A12u[1] = A12u[1] + pA12[i*3+1]*u[i];
        A12u[2] = A12u[2] + pA12[i*3+2]*u[i];
        
        // A12*v [3x3] 
	    A12v[0] = A12v[0] + pA12[i*3+0]*v[i*3+0] ;
        A12v[1] = A12v[1] + pA12[i*3+0]*v[i*3+1] ;      
        A12v[2] = A12v[2] + pA12[i*3+0]*v[i*3+2] ;       
        
        A12v[1*3+0] = A12v[1*3+0] + pA12[i*3+1]*v[i*3+0] ;
        A12v[1*3+1] = A12v[1*3+1] + pA12[i*3+1]*v[i*3+1] ;
        A12v[1*3+2] = A12v[1*3+2] + pA12[i*3+1]*v[i*3+2] ;
        
        A12v[2*3+0] = A12v[2*3+0] + pA12[i*3+2]*v[i*3+0] ;
        A12v[2*3+1] =  A12v[2*3+1] + pA12[i*3+2]*v[i*3+1] ;
        A12v[2*3+2] = A12v[2*3+2] + pA12[i*3+2]*v[i*3+2] ;    
          
	}

     
    A_[0] = pA11[0] - A12v[0];
    A_[1] = pA11[1] - A12v[1];
    A_[2] = pA11[2] - A12v[2];
    
 
    
    A_[3+0] = pA11[3+0] - A12v[3+0];
    A_[3+1] = pA11[3+1] - A12v[3+1];
    A_[3+2] = pA11[3+2] - A12v[3+2];
    
    A_[2*3+0] = pA11[2*3+0] - A12v[2*3+0];
    A_[2*3+1] = pA11[2*3+1] - A12v[2*3+1];
    A_[2*3+2] = pA11[2*3+2] - A12v[2*3+2];
           
    // b(1:3) - A12u
    b_[0] = b[0] - A12u[0];
    b_[1] = b[1] - A12u[1];
    b_[2] = b[2] - A12u[2];
    

    inverse3x3Matrix(A_inv, A_);

    
   x[0] = A_inv[0*3+0]*b_[0] + A_inv[0*3+1]*b_[1] + A_inv[0*3+2]*b_[2] ;
   x[1] = A_inv[1*3+0]*b_[0] + A_inv[1*3+1]*b_[1] + A_inv[1*3+2]*b_[2] ;
   x[2] = A_inv[2*3+0]*b_[0] + A_inv[2*3+1]*b_[1] + A_inv[2*3+2]*b_[2] ;
    
   
   for (int i=0;i<n2;i++)
   {
       double temp = (v[i*3+0]*x[0] + v[i*3+1]*x[1] + v[i*3+2]*x[2]);
    
        x[i+3] = u[i] - temp;
   }
   
   delete [ ] u;
   delete [ ] v;

	return;
}

void inverse3x3Matrix(double *A_i, const double *A_)
{
    if (A_ == 0 || A_i == 0)
        return;
       // 3x3 inverse 
    double A_det = (A_[0]*A_[1*3+1]*A_[2*3+2] - A_[0]*A_[1*3+2]*A_[2*3+1] - A_[1]*A_[1*3]*A_[2*3+2] + A_[1]*A_[1*3+2]*A_[2*3] + A_[2]*A_[1*3]*A_[2*3+1] - A_[2]*A_[1*3+1]*A_[2*3]);
    
    A_i[0] =   (A_[1*3+1]*A_[2*3+2] - A_[1*3+2]*A_[2*3+1])/A_det;
    A_i[1] = -(A_[0*3+1]*A_[2*3+2] - A_[0*3+2]*A_[2*3+1])/A_det;
    A_i[2] =   (A_[0*3+1]*A_[1*3+2] - A_[0*3+2]*A_[1*3+1])/A_det;
    
    A_i[3+0] =  -(A_[1*3+0]*A_[2*3+2] - A_[1*3+2]*A_[2*3+0])/A_det;
    A_i[3+1] =   (A_[0*3+0]*A_[2*3+2] - A_[0*3+2]*A_[2*3+0])/A_det;
    A_i[3+2] =  -(A_[0*3+0]*A_[1*3+2] - A_[0*3+2]*A_[1*3+0])/A_det;
    
    A_i[2*3+0] =   (A_[1*3+0]*A_[2*3+1] - A_[1*3+1]*A_[2*3+0])/A_det;
    A_i[2*3+1] =  -(A_[0*3+0]*A_[2*3+1] - A_[0*3+1]*A_[2*3+0])/A_det;
    A_i[2*3+2] =   (A_[0*3+0]*A_[1*3+1] - A_[0*3+1]*A_[1*3+0])/A_det;    
        
}

