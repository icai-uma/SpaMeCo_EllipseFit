// Sample some points from an ellipse in parametric form
// Original C++ version:
// void CalculateEllipsePoints(ParamEllipse pe, Pixel* points, int noPoints);
// MEX version:
// [points]=CalculateEllipsePointsMEX(MyEllipse,noPoints);
// where MyEllipse=[xC, yC, a, b, rotationDegree]
// MEX file to be compiled with this command at the Matlab prompt:
// mex CalculateEllipsePointsMEX.cpp EECM.cpp
// Conversion to Matlab MEX by Ezequiel Lopez-Rubio
// Version Date: 14th August 2020

#include <iostream>
#include "mex.h"
#include "EECM.h"

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* MyEllipse;
    int noPoints;
    double *ptrOutput;
    
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "CalculateEllipsePointsMEX requires two input arguments.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "CalculateEllipsePointsMEX requires one output argument.");
    }
    MyEllipse = mxGetPr(prhs[0]);
    noPoints = mxGetScalar(prhs[1]);
    
    Pixel* points = new Pixel[noPoints];
    ParamEllipse pe(MyEllipse[0],MyEllipse[1],MyEllipse[2],MyEllipse[3],MyEllipse[4]);
    CalculateEllipsePoints(pe, points, noPoints);
    plhs[0]=mxCreateDoubleMatrix(2, noPoints, mxREAL);
    ptrOutput=mxGetPr(plhs[0]);
    for(int NdxPoint=0;NdxPoint<noPoints;NdxPoint++)
    {
        ptrOutput[2*NdxPoint]=points[NdxPoint].x;
        ptrOutput[2*NdxPoint+1]=points[NdxPoint].y;
    }
    delete points;
}
