#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "mex.h"
#include "matrix.h"
/*******************************************************************/
/* 

 fast sigmoid. Used in sigmoid_approx
 
%%%%% BEWARE %%%%%
MATLAB code to compile it : mex mex_sigmoid_fast.c

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

*/

        
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* input */    
    double *x_l;
    x_l = mxGetData(prhs[0]);

    int N;
    N = mxGetNumberOfElements(prhs[0]);
    
    /* output */
    double *y;
    plhs[0]=mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
    y = mxGetData(plhs[0]);

    /*  additional parameters*/
    
    double xbnd_l[] = {-2.0, -0.5, 0.5, 2}; /* boundaries  */
    int fbnd = 3; /* index of final bound (length(xbnd_l) -1)*/
    double x0_l[] = {-1, 0, 1}; /* reference values */
    double d0_l[] = {0.2689,  0.5, 0.7311 }; /* values */
    double d1_l[] = {0.1966, 0.25, 0.1966}; /* 1st derivative */
    double d2_l[] = {0.0909/2.0, 0.0, -0.0909/2.0}; /* 2nd derivative  /2 */
    
    /* computation */
    int i, f;
    double x, d;
    for (i=0;i<N;i++){
        x = x_l[i];
        if (x< xbnd_l[0]){
            y[i] = 1/(1+exp(-x));
        }
        else if (x<xbnd_l[fbnd]){
            f = fbnd-1;
            while (x < xbnd_l[f+1]){
                f = f-1;                
            }
            f = f+1;
            d = x - x0_l[f];
            y[i] = d0_l[f] + d*d1_l[f] + d*d*d2_l[f];
        }
        else {
            y[i] = 1/(1+exp(-x));
        }       
    }        
	return;		
}
