// This file declares Euclidean Ellipse Comparison Metric (EECM),
// an ellipse fit error and ellipse distance calculation method based on
// the method in [1]. Cite this work as:
// 
// Cakir, H. I., Topal, C. "An Euclidean Ellipse Comparison Metric 
// for Quantitative Evaluation", International Conference on 
// Acoustics, Speech and Signal Processing (ICASSP), 2018. 
// 
// Note that this file is designed to work with the "Screen Coordinate 
// System" which has the origin at the top-left corner.
// This is why the word "Pixel" is especially preferred instead of 
// "Point". Thus, all points' y components are reflected with a 
// simple transform.
// 
// [1] R. Nurnberg, Distance from a Point to an Ellipse, 2006.
// Online available at www.ma.ic.ac.uk/~rn/distance2ellipse.pdf.
//
// Version Date: 2018-02-15
// Halil Ibrahim Cakir & Cihan Topal
// For more information: http://c-viz.anadolu.edu.tr/eecm

#ifndef EECM_H
#define EECM_H

#define _USE_MATH_DEFINES

#include <math.h> // defines M_PI

struct Pixel
{

public:

	double x, y;

	Pixel(){}

	Pixel(double x, double y)
	{
		this->x = x;
		this->y = y;
	}
};

// Conic Equation of Ellipse
struct ConicEllipse
{

public:

	double A, B, C, D, E, F;

	ConicEllipse(double A, double B, double C, double D, double E, double F)
	{
		this->A = A;
		this->B = B;
		this->C = C;
		this->D = D;
		this->E = E;
		this->F = F;
	}
};

// Parametric Equation of Ellipse
struct ParamEllipse
{

private:

	double xC, yC, a, b, rotationDegree, rotationRadian;

	double a2_b2;

public:

	ParamEllipse(double xC, double yC, double a, double b, double rotationDegree)
	{
		this->xC = xC;
		this->yC = yC;
		this->a = a;
		this->b = b;
		this->rotationDegree = rotationDegree;
		this->rotationRadian = rotationDegree * M_PI / 180.0;

		this->a2_b2 = a * a - b * b;
	}

	double GetCenterX(){ return xC; }
	double GetCenterY(){ return yC; }
	double GetSemiMajor(){ return a; }
	double GetSemiMinor(){ return b; }
	double GetRotationDegree(){ return rotationDegree; }
	double GetRotationRadian(){ return rotationRadian; }

	double Geta2_b2(){ return a2_b2; }
};

void CalculateEllipsePoints(ParamEllipse pe, Pixel* points, int noPoints);

double NewtonRaphsonIteration(double theta, double x, double y, ParamEllipse pe);

double NewtonRaphsonThetaEstimation(double th0, double x, double y, ParamEllipse e);

double DistanceToPoint(ParamEllipse pe, Pixel point);

double DistanceBwEllipses(ParamEllipse pe1, ParamEllipse pe2, int noPoints);

ParamEllipse ConvertEllipse(ConicEllipse ce);

#endif