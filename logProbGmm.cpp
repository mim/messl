/*
 * lp = logProbGmm(Erep, xi_wit, s2_wit, log(s2_wit));
 *
 * Compute the log probability of Erep being in a Gaussian with means
 * xi_wit and variances s2_wit.  Erep is a 4D array with size [W T I
 * Nt], xi_wit and s2_wit are both 3D arrays with size [W I Nt] that
 * need to be replicated T times along Erep's second dimension.  The
 * fourth argument is the log of the third because I was having
 * trouble with memory allocation if I computed it inside this
 * function.
 */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ 
  if(nrhs != 4)
    mexErrMsgTxt("Usage: lp = logProbGmm(Erep, xi_wit, s2_wit)");
  if(!mxIsSingle(prhs[0]))
    mexErrMsgTxt("First arguments must be single");
  if(!(mxIsDouble(prhs[1]) && mxIsDouble(prhs[2]) && mxIsDouble(prhs[3])))
    mexErrMsgTxt("All arguments after the first must be doubles");

  if(mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) || 
     mxIsComplex(prhs[2]) || mxIsComplex(prhs[3]))
    mexErrMsgTxt("Arguments must be real");

  if(mxGetNumberOfDimensions(prhs[0]) != 4)
     mexErrMsgTxt("Erep must have 4 dimensions");
  if(mxGetNumberOfDimensions(prhs[1]) != 3)
    mexErrMsgTxt("xi_wit must have 3 dimensions");
  if(mxGetNumberOfDimensions(prhs[2]) != 3)
    mexErrMsgTxt("s2_wit must have 3 dimensions");
  if(mxGetNumberOfDimensions(prhs[3]) != 3)
    mexErrMsgTxt("log(s2_wit) must have 3 dimensions");

  /* Should be mwSize type? */
  const int *dims = mxGetDimensions(prhs[0]);
  int W=dims[0], T=dims[1], I=dims[2], Nt=dims[3];

  /* TODO: test that mu_wit and s2_wit are the right size */

  float  *Erep   = (float*)mxGetData(prhs[0]);
  double *xi_wit = mxGetPr(prhs[1]);
  double *s2_wit = mxGetPr(prhs[2]);
  double *logS2  = mxGetPr(prhs[3]);

  /* Easiest way to make a 4D matrix of the same size */
  plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
  double *lp = mxGetPr(plhs[0]);

  int w,t,i,tau;    /* Should be mwIndex type? */

  /* TODO: Make an array that holds the log of s2_wit. */
/*   double *logS2 = (double *)mxCalloc(W*I*Nt, sizeof(double)); */
/*   for (i = 0; i < W*I*Nt; i++) */
/*     logS2[i] = log(s2_wit[i]); */

  for (tau = 0; tau < Nt; tau++) {
    for (i = 0; i < I; i++) {
      for (t = 0; t < T; t++) {
        for (w = 0; w < W; w++) {
          int i4d = w + W*t + W*T*i + W*T*I*tau;
          int i3d = w + W*i + W*I*tau;
          double diff = Erep[i4d] - xi_wit[i3d];
          lp[i4d] = -0.5*logS2[i3d] - 0.5 *diff*diff / s2_wit[i3d];
        }
      }
    }
  }
}
