function [ParN]=GtoN(ParG)
% Compute the natural parameters of an ellipse from the geometric parameters
% Input:
%   ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]' is the vector of 
%   geometric parameters of the ellipse). a=half length of the major axis,
%   b=half length of the minor axis
% Output:
%   ParN= [Focus1x Focus1y Focus2x Focus2y SumDists]' is the vector of
%   natural parameters of the ellipse: the two foci (Focus1 and Focus2),
%   and the sum of distances to both foci, SumDists==2a

Center=ParG(1:2);
a=ParG(3);
b=ParG(4);
phi=-ParG(5);

CenterToFocusDistance=sqrt(a^2-b^2);
Focus1=[-CenterToFocusDistance 0]';
Focus2=[CenterToFocusDistance 0]';
RotationMatrix=[cos(-phi) -sin(-phi);sin(-phi) cos(-phi)];
Focus1=RotationMatrix*Focus1+Center;
Focus2=RotationMatrix*Focus2+Center;

ParN=[Focus1(1) Focus1(2) Focus2(1) Focus2(2) 2*a]';
