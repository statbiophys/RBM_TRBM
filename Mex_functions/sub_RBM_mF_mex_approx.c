#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"
/*******************************************************************/
/*
%%%%% BEWARE %%%%%
MATLAB code to compile it : mex sub_RBM_mF_mex_approx.c

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

*/

        
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* input */
    
    double *a, *b, *w;
    a=mxGetData(mxGetField(prhs[0], 0, "a"));
    b=mxGetData(mxGetField(prhs[0], 0, "b"));
    w=mxGetData(mxGetField(prhs[0], 0, "w")); 

    double *wv;
    wv = mxGetData(prhs[1]);
    
    double *x_max_D, x_max;
    x_max_D = mxGetData(prhs[2]);
    x_max = x_max_D[0];
    
    int Ni, Nj, Nt;
    Nj = mxGetM(mxGetField(prhs[0], 0, "w"));
    Ni = mxGetN(mxGetField(prhs[0], 0, "w"));
    Nt = mxGetN(prhs[1]);
    
     /* check input  */
    if (Nj != mxGetM(prhs[1])) { 
        mexErrMsgTxt("M.w*v and M.w have different Nj: check second argument is M.w*v and not v alone\n");
    }
          
    /* output */
    double *mF;
 
    plhs[0]=mxCreateDoubleMatrix(1,Nt,mxREAL);
    mF = mxGetData(plhs[0]);

    int t, j;
    double x, pro;
    for (t=0;t<Nt;t++){
        pro = 1;
        for (j=0;j<Nj;j++){
            x = wv[j+t*Nj]+b[j];
            if (x>x_max) {
                mF[t] = mF[t] +x;
            }
            else {
                pro = pro*(1+exp(x));
            }
        }
        mF[t] = mF[t]+log(pro);
    }
          
	return;		
}
