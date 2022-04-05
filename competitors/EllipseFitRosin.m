function [ParA,ParG,ParN,code]=EllipseFitRosin(X)
% Rosin method, normalization A+C=1
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


x=X(1,:);
y=X(2,:);
B = [sum(((x.^2)-(y.^2)).^2) sum(((x.^2)-(y.^2)).*x.*y) sum(((x.^2)-(y.^2)).*(x)) sum(((x.^2)-(y.^2)).*(y)) sum(((x.^2)-(y.^2))); ...
    sum(x.*y.*((x.^2)-(y.^2))) sum((x.^2).*(y.^2)) sum((x.^2).*(y)) sum((x).*(y.^2)) sum((x).*(y)); ...
    sum(x.*((x.^2)-(y.^2))) sum((x.^2).*(y)) sum(x.^2) sum((x).*(y)) sum(x); ...
    sum(y.*((x.^2)-(y.^2))) sum((x).*(y.^2)) sum((x).*(y)) sum((y.^2)) sum((y)); ...
    sum(((x.^2)-(y.^2))) sum((x).*(y)) sum(x) sum(y) sum(ones(size(x))); ...
    ];
C = [-(sum(((x.^2)-(y.^2)).*(y.^2))) -(sum(x.*(y.^3))) -(sum(x.*(y.^2))) -(sum(y.^3)) -(sum(y.^2))];
A =(C*inv(B));
%a2=A(1);b2=A(2);c2=1-A(1);d2=A(3);e2=A(4);f2=A(5);

ParA(1) = A(1);
ParA(2) = 0.5*A(2);
ParA(3) = 1-A(1);
ParA(4) = 0.5*A(3);
ParA(5) = 0.5*A(4);
ParA(6) = A(5);

ParA = ParA';

[ParG,code]=AtoG(ParA);
if code==8, ParN=NaN(5,1);  % If the fit is a single point GtoN fails
else,ParN=GtoN(ParG);
end
