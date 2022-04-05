function [ParA,ParG,ParN]=EllipseFitSzpak(X)
% Szpak method
% Inputs:
%   X=matrix of size 2 x NumPoints with the training samples
% Outputs:
%   ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]' is the vector of 
%   geometric parameters of the ellipse). a=half length of the major axis,
%   b=half length of the minor axis
%   ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters:
%           Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0
%   ParN= [Focus1x Focus1y Focus2x Focus2y SumDists] is the vector of
%   natural parameters of the ellipse: the two foci (Focus1 and Focus2),
%   and the sum of distances to both foci, SumDists==2a


ParA = fast_guaranteed_ellipse_estimate(X');

% geometricEllipseParameters = ...
%             fromAlgebraicToGeometricParameters(ParA);
% ParG=[geometricEllipseParameters(3:4)' geometricEllipseParameters(1:2)' geometricEllipseParameters(5)]';

ParG=AtoG([ParA(1) ParA(2)/2 ParA(3) ParA(4)/2 ParA(5)/2 ParA(6)]);
ParN=GtoN(ParG);
