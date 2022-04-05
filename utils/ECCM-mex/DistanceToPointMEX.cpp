// Compute distance from an ellipse in parametric form to a point
// Original C++ version:
// double DistanceToPoint(ParamEllipse pe, Pixel point);
// MEX version:
// MyDistance=DistanceToPointMEX(MyEllipse, MyPixel);
// where MyEllipse=[xC, yC, a, b, rotationDegree], MyPixel=[x,y]
// MEX file to be compiled with this command at the Matlab prompt:
// mex DistanceToPointMEX.cpp EECM.cpp
// Conversion to Matlab MEX by Ezequiel Lopez-Rubio
// Version Date: 14th August 2020

#include <iostream>
#include "mex.h"
#include "EECM.h"

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* MyEllipse;
    double* MyPoint;
    double *ptrOutput;
    
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "DistanceToPointMEX requires two input arguments.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "DistanceToPointMEX requires one output argument.");
    }
    MyEllipse = mxGetPr(prhs[0]);
    MyPoint = mxGetPr(prhs[1]);
    
    ParamEllipse pe(MyEllipse[0],MyEllipse[1],MyEllipse[2],MyEllipse[3],MyEllipse[4]);
    Pixel point(MyPoint[0],MyPoint[1]);   
    plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
    ptrOutput=mxGetPr(plhs[0]);
    
    (*ptrOutput)=DistanceToPoint(pe, point);

}
