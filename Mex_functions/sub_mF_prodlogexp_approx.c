#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"
/*******************************************************************/
/*
MATLAB code to compile it : mex sub_mF_prodlogexp_approx.c

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017
 
*/

        
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* input */
    double *wv;
    wv = mxGetData(prhs[0]);
    
    double *x_max_D, x_max;
    x_max_D = mxGetData(prhs[1]);
    x_max = x_max_D[0];
    
    int Nj, Nt;
    Nj = mxGetM(prhs[0]);
    Nt = mxGetN(prhs[0]);
          
    /* output */
    double *mF;
 
    plhs[0]=mxCreateDoubleMatrix(1,Nt,mxREAL);
    mF = mxGetData(plhs[0]);

    int t, j;
    double x, pro;
    for (t=0;t<Nt;t++){
        pro = 1;
        for (j=0;j<Nj;j++){
            x = wv[j+t*Nj];
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
