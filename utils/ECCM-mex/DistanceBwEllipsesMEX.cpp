// Compute distance from an ellipse in parametric form to a point
// Original C++ version:
// double DistanceBwEllipses(ParamEllipse pe1, ParamEllipse pe2, int noPoints);
// MEX version:
// MyDistance=DistanceBwEllipsesMEX(MyEllipse1, MyEllipse2, noPoints);
// where MyEllipse1=[xC1, yC1, a1, b1, rotationDegree1], MyEllipse2=[xC2, yC2, a2, b2, rotationDegree2]
// MEX file to be compiled with this command at the Matlab prompt:
// mex DistanceBwEllipsesMEX.cpp EECM.cpp
// Conversion to Matlab MEX by Ezequiel Lopez-Rubio
// Version Date: 14th August 2020

#include <iostream>
#include "mex.h"
#include "EECM.h"

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* MyEllipse1;
    double* MyEllipse2;
    int noPoints;
    double *ptrOutput;
    
    
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "DistanceBwEllipsesMEX requires two input arguments.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "DistanceBwEllipsesMEX requires one output argument.");
    }
    MyEllipse1 = mxGetPr(prhs[0]);
    MyEllipse2 = mxGetPr(prhs[1]);
    noPoints = mxGetScalar(prhs[2]);
    
    ParamEllipse pe1(MyEllipse1[0],MyEllipse1[1],MyEllipse1[2],MyEllipse1[3],MyEllipse1[4]);
    ParamEllipse pe2(MyEllipse2[0],MyEllipse2[1],MyEllipse2[2],MyEllipse2[3],MyEllipse2[4]);  
    plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
    ptrOutput=mxGetPr(plhs[0]);
    
    (*ptrOutput)=DistanceBwEllipses(pe1, pe2, noPoints);

}
