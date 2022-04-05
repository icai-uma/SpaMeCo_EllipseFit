// This file defines Euclidean Ellipse Comparison Metric (EECM),
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

// Conversion to Matlab MEX by Ezequiel Lopez-Rubio
// Version Date: 14th August 2020

#include "EECM.h"

// Calculates sampled points on ellipse contour
void CalculateEllipsePoints(ParamEllipse pe, Pixel* points, int noPoints)
{
	double xC = pe.GetCenterX(); double yC = pe.GetCenterY(); double a = pe.GetSemiMajor(); double b = pe.GetSemiMinor();

	// in Cartesian
	double tmpYc = -yC; // y-axis is reflected from + to -
	double tmpRotationRadian = -pe.GetRotationRadian(); // rotation converted from clockwise to counter-clockwise

	float angleStep = 360.0 / (float)noPoints;
	float currentAngle = 0.0;
	double tmpX, tmpY;

	for (int i = 0; i < noPoints; i++)
	{
		// assume that ellipse is in origin and unrotated
		points[i].x = a * sin(currentAngle);
		points[i].y = b * cos(currentAngle);

		// keep point temporarily
		tmpX = points[i].x;
		tmpY = points[i].y;

		// rotate point
		points[i].x = tmpX * cos(tmpRotationRadian) - tmpY * sin(tmpRotationRadian);
		points[i].y = tmpX * sin(tmpRotationRadian) + tmpY * cos(tmpRotationRadian);

		// translate point to origin
		points[i].x += xC;
		points[i].y += tmpYc; // in Cartesian

		points[i].y *= -1; // in Screen

		currentAngle += angleStep;
	}
}

// Returns (F / F') for one NR iteration
double NewtonRaphsonIteration(double theta, double x, double y, ParamEllipse pe)
{
	double a = pe.GetSemiMajor(), b = pe.GetSemiMinor(), a2_b2 = pe.Geta2_b2();

	// compute the frequently used terms for the sake of performance
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);

	double fTheta = a2_b2 * cosTheta * sinTheta;
	fTheta = fTheta - x * a * sinTheta + y * b * cosTheta;

	double fPrimeTheta = a2_b2 * (cosTheta*cosTheta - sinTheta*sinTheta);
	fPrimeTheta = fPrimeTheta - x * a * cosTheta - y * b * sinTheta;

	return (fTheta / fPrimeTheta);
}

// The number of iterations can be varied for higher accuracies
double NewtonRaphsonThetaEstimation(double th0, double x, double y, ParamEllipse e)
{
	double th1 = th0 - NewtonRaphsonIteration(th0, x, y, e);
	double th2 = th1 - NewtonRaphsonIteration(th1, x, y, e);
    //Nuevas
    double th3 = th2 - NewtonRaphsonIteration(th2, x, y, e);
    double th4 = th3 - NewtonRaphsonIteration(th3, x, y, e);

	//return th2;
    return th4;
}

// Calculates distance from point to ellipse
double DistanceToPoint(ParamEllipse pe, Pixel point)
{
	double xC = pe.GetCenterX(), yC = pe.GetCenterY(), a = pe.GetSemiMajor(), b = pe.GetSemiMinor(),
		rotationRadian = pe.GetRotationRadian();

	// different initial values for 4 estimations (in radians)
	double th0[4], estimations[4];

	// closest point candidates and distances for different estimations
	double distances[4], x, y;
	double candidatesX[4], candidatesY[4];
	double transPntX, transPntY;	// translated and rotated Point

	// translate & rotate candidate point coordinates
	transPntX = point.x - xC;
	transPntY = point.y - yC;

	x = transPntX * cos(-rotationRadian) - transPntY * sin(-rotationRadian); // rotation is in reverse direction
	y = transPntX * sin(-rotationRadian) + transPntY * cos(-rotationRadian);

	transPntX = x;
	transPntY = y;

	th0[0] = atan(a*x / b*y);
    th0[1] = th0[0] + M_PI*0.5;	//pi * 0.5;
	th0[2] = th0[0] + M_PI;		//pi * 1.0
	th0[3] = th0[0] + M_PI*1.5;	//pi * 1.5
	/*
    th0[0] = atan(a*x / a*y);
    th0[1] = th0[0] + 1.5707;	//pi * 0.5;
	th0[2] = th0[0] + M_PI;		//pi * 1.0
	th0[3] = th0[0] + 4.7123;	//pi * 1.5
     */

	int minDistPnt = 0;

	for (int i = 0; i < 4; i++)
	{
		estimations[i] = NewtonRaphsonThetaEstimation(th0[i], x, y, pe);

		candidatesX[i] = a * cos(estimations[i]);
		candidatesY[i] = b * sin(estimations[i]);

		distances[i] = sqrt((transPntX - candidatesX[i])*(transPntX - candidatesX[i]) + (transPntY - candidatesY[i])*(transPntY - candidatesY[i]));

		if (distances[i] < distances[minDistPnt]) minDistPnt = i;
	}

	double tmpX, tmpY;

	tmpX = a * cos(estimations[minDistPnt]);
	tmpY = b * sin(estimations[minDistPnt]);

	return distances[minDistPnt];
}

// Calculates distance between two ellipses
double DistanceBwEllipses(ParamEllipse pe1, ParamEllipse pe2, int noPoints)
{
	/* Calculate ellipse drawing points */

	Pixel* pe1Points = new Pixel[noPoints];
	Pixel* pe2Points = new Pixel[noPoints];

	CalculateEllipsePoints(pe1, pe1Points, noPoints);
	CalculateEllipsePoints(pe2, pe2Points, noPoints);

	/* Calculate distance from each ellipse point to other ellipse */

	double sumDistance = 0.0;

	for (int i = 0; i < noPoints; i++)
	{
		sumDistance += DistanceToPoint(pe2, pe1Points[i]) + DistanceToPoint(pe1, pe2Points[i]);
	}

	delete[] pe1Points;
	delete[] pe2Points;

	return sumDistance / (double)(2 * noPoints);
}

// Converts conic ellipse to parametric one
ParamEllipse ConvertEllipse(ConicEllipse ce)
{
	double xC, yC, a, b, rotation;

	double A1 = ce.A, A2;
	double B1 = ce.B, B2;
	double C1 = ce.C, C2;
	double D1 = ce.D, D2;
	double E1 = ce.E, E2;
	double F1 = ce.F, F2;

	/* Convert the conic equation to the ellipse equation */
	
	// normalize coefficients
	B1 /= A1; C1 /= A1; D1 /= A1; E1 /= A1; F1 /= A1; A1 /= A1;

	if (B1 == 0) // then not need to rotate the axes
	{
		A2 = A1;
		B2 = B1;
		C2 = C1;
		D2 = D1;
		E2 = E1;
		F2 = F1;
	}
	else if (B1 != 0) // rotate the axes
	{
		/* Determine the rotation angle (in radians) */

		if (B1 == 0 && A1 < C1)
			rotation = 0;
		else if (B1 == 0 && A1 > C1)
			rotation = M_PI / 2;
		else if (B1 != 0 && A1 < C1)
			rotation = atan(B1 / (A1 - C1)) / 2;
		else if (B1 != 0 && A1 > C1)
			rotation = M_PI / 2 + atan(B1 / (A1 - C1)) / 2;
		else if (B1 != 0 && A1 == C1)
			rotation = M_PI / 4;
		else
			rotation = atan(B1 / (A1 - C1)) / 2;

		// compute the coefficients wrt the new coordinate system 
		A2 = 0.5 * (A1 * (1 + cos(2 * rotation) + B1 * sin(2 * rotation) + C1 * (1 - cos(2 * rotation))));
		B2 = (C1 - A1) * sin(2 * rotation) + B1 * cos(2 * rotation); //B2 should turn to be zero?
		C2 = 0.5 * (A1 * (1 - cos(2 * rotation) - B1 * sin(2 * rotation) + C1 * (1 + cos(2 * rotation))));
		D2 = D1 * cos(rotation) + E1 * sin(rotation);
		E2 = -D1 * sin(rotation) + E1 * cos(rotation);
		F2 = F1;
	}

	// transform the conic equation into the ellipse form
	double D3 = D2 / A2; // normalize x term's coef 

	double E3 = E2 / C2; // normalize y term's coef

	xC = -(D3 / 2);	// center X
	yC = -(E3 / 2);	// center Y

	double F3 = A2 * (xC * xC) + C2 * (yC * yC) - F2;

	// semimajor axis
	a = sqrt(F3 / A2);
	// semiminor axis
	b = sqrt(F3 / C2);

	// center coordinates have to be re-transformed if rotation is applied!
	if (rotation != 0)
	{
		double tmpX = xC, tmpY = yC;
		xC = tmpX * cos(rotation) - tmpY * sin(rotation);
		yC = tmpX * sin(rotation) + tmpY * cos(rotation);
	}

	return ParamEllipse(xC, -yC, a, b, -rotation * 180.0 / M_PI);
}