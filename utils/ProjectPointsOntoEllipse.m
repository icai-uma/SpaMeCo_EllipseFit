function [XYproj,ParameterProj] = ProjectPointsOntoEllipse(XY,ParG)
% http://people.cas.uab.edu/~mosya/cl/MATLABconics.html
%   Projecting a given set of points onto an ellipse
%   and computing the distances from the points to the ellipse
%
%   This is a modified version of an iterative algorithm published by D. Eberly
%     Internet publication: "Distance from a point to an ellipse in 2D" (2004)
%                           Geometric Tools, LLC, www.geometrictools.com
%     Book publication: "3D Game Engine Design", 2nd edition.
%                       Morgan Kaufmann Publishers, San Francisco, CA, 2007.
%                              (see Section 14.13.1)
%
%   Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%           ParG is a vector 5x1 of the ellipse parameters
%           ParG =  [Center(1:2), Axes(1:2), Angle]'
%                Center - the coordinates of the ellipse's center
%                Axes   - the axes (major, minor)
%                Angle  - the angle of tilt of the ellipse
%
%           Note: the user needs to make sure that ParG(2) >= ParG(3) > 0.
%
%   Output:  XYproj is the array of coordinates of projections
%           ParameterProj is the array of parameters t of the projections, so
%           that in the canonical coordinate system we have:
%           XYproj=(a*cos(t),b*sin(t))
%
%   The algorithm is proven to converge and reaches an accuracy of 10-12 significant digit
%   It takes 4-6 iterations per point, on average.
%
%        Nikolai Chernov,  March 2012
% Note: in the canonical coordinate system, the ellipse passes through the
% points (a,0) and (0,b)

Center = ParG(1:2);   a = ParG(3);  b = ParG(4);   Angle = ParG(5);  % ellipse parameters

if (a <= 0)||(b <= 0)
    disp(' Axes of the ellipse must be positive')
    XYproj = NaN;
    ParameterProj = NaN;
    return
end
if a < b
    disp(' Major axis of the ellipse cannot be smaller than its minor axis')
    XYproj = NaN;
    ParameterProj = NaN;
    return
end

n = size(XY,1);       %  n is the number of points to be projected
XYproj = zeros(n,2);  %  prepare an array of projections

aa = a^2;  bb = b^2;          %  squares of the axes
D  = (a-b)*(a+b);             %  "distorsion measure"

Q = [cos(Angle) -sin(Angle); sin(Angle) cos(Angle)];
%    Q is the matrix for rotating the points and the hyperbola to the canonical system

XY0  = [XY(:,1)-Center(1) XY(:,2)-Center(2)]*Q;    %  points in canonical coordinates

for i=1:n        %  Main loop over the points to be projected
    
    u = abs(XY0(i,1));  v = abs(XY0(i,2));     %  coordinates of the point
    T = max(a*u-D,b*v);        %  initial value of the T variable
    if (T <= 0)&&(D <= 0)    %  circle (a=b) and point at its center
        XYproj(i,1) = 0;  XYproj(i,2) = b;
        continue;
    end
    if T <= 0      %  true ellipse (a>b) and point on major axis near or at center
        XYproj(i,1) = aa*XY0(i,1)/D;  XYproj(i,2) = b*sqrt(max(1-(a*u/D)^2,0));
        continue;
    end
    
    %  the main, non-singular case
    %  start Newton's iterations.
    
    iterMax = 100; %   Maximum number of Newton's ietrations. Usually, 4-6 are enough
     
    for iter=1:iterMax %  loop of Newton's iterations to solve F(T)=0
        F  = (a*u/(T+D))^2 + (b*v/T)^2 - 1;      %  value of F; we need to find T such that F(T)=0
        if F <= 0                %   gone too far, emergency stop
            break;
        end;
        Fder = -2*((a*u/(T+D))^2/(T+D) + (b*v/T)^2/T);  %  derivative of F with respect to T
        Step = F/Fder;
        if T == T - Step, break; end;   %  no progress, terminate iterations
        T = T - Step;                   %  Newton's iteration
    end    %  end Newton's iterations
    
    %            compute the projection of the point onto the ellipse
    
    xprojx = aa*u/(T+D);            %   first candidate for projection
    yprojx = b*sqrt(max(1-(xprojx/a)^2,0));
    yprojy = bb*v/T;                %   second candidate for projection
    xprojy = a*sqrt(max(1-(yprojy/b)^2,0));
    
    Fx = (xprojx-u)^2 + (yprojx-v)^2;  % "quality" of first candidate
    Fy = (xprojy-u)^2 + (yprojy-v)^2;  % "quality" of second candidate
    
    if Fx < Fy         %    the first  candidate is better
        XYproj(i,1) = xprojx;
        XYproj(i,2) = yprojx;
    else               %    the second candidate is better
        XYproj(i,1) = xprojy;
        XYproj(i,2) = yprojy;
    end                %    end comparing the two candidates
    
    %         compute the projection point in the proper quadrant
    
    if (XY0(i,1)<0),  XYproj(i,1) = -XYproj(i,1);  end;
    if (XY0(i,2)<0),  XYproj(i,2) = -XYproj(i,2);  end;

end     %  end the main loop over the points

ParameterProj=atan2(XYproj(:,2)/b,XYproj(:,1)/a);

%    rotate and shift back to the original system

XYproj = XYproj*Q';     %  rotation
XYproj(:,1) = XYproj(:,1)+Center(1);  %  shifting
XYproj(:,2) = XYproj(:,2)+Center(2);  %  shifting

XYproj=XYproj';
end   %  ProjectPointsOntoEllipse
