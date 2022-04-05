function ParA = GtoA(ParG,code)
% http://people.cas.uab.edu/~mosya/cl/MATLABconics.html

%  Conversion of Geometric parameters of a conic to its Algebraic parameters

%   Geometric parameters depend on the type of the conic (ellipse, etc.)

%   Algebraic parameters are coefficients A,B,C,D,E,F in the algebraic
%   equation     Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0

%   Input:  code is the code of the conic type (see below)
%           ParG is the vector of Geometric parameters (see below)

%   Output: ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters

%   code:   1 - ellipse  2 - hyperbola  3 - parabola
%           4 - intersecting lines  5 - parallel lines
%           6 - coincident lines    7 -single line
%           8 - single point        9 - imaginary ellipse
%          10 - imaginary parallel lines
%          11 - "impossible" equation, 1=0 or -1=0 (no solutions)

%           Geometric parameters are determined only for 
%               the most common types of conics:

%   1. Ellipses:  canonical equation  x^2/a^2 + y^2/b^2 = 1
%                 a>b or a=b are the major and minor semi-axes

%       ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]'

%   2. Hyperbolas:  canonical equation  x^2/a^2 - y^2/b^2 = 1
%                 a  is the distance from the center to each vertex
%                 b  is the vertical distance from each vertex to asymptotes

%       ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]'

%   3. Parabolas:  canonical equation  y^2 = 2px
%                 p  is the distance from the focus to the directix

%       ParG = [Xcenter, Ycenter, p, AngleOfTilt]'

%   4. Intersecting lines:  canonical equation  
%                (cos(u)*x+sin(u)*y+d)*(cos(v)*x+sin(v)*y+e) = 0
%                   u,v are directional angles of the lines
%                   d,e are distances from the origin (0,0) to the lines

%       ParG = [u,d,v,e]'

%        Nikolai Chernov,  February 2012

if (code == 1)                           %   ellipse
    c = cos(ParG(5));  s = sin(ParG(5));
    a = ParG(3);  b = ParG(4);
    Xc = ParG(1);  Yc = ParG(2);
    P = (c/a)^2 + (s/b)^2;
    Q = (s/a)^2 + (c/b)^2;
    R = c*s*(1/a^2 - 1/b^2);
    ParA = [P; R; Q; -P*Xc-R*Yc; -Q*Yc-R*Xc; P*Xc^2+Q*Yc^2+2*R*Xc*Yc-1];
    ParA = ParA/norm(ParA);
end

if (code == 2)                           %   hyperbola
    c = cos(ParG(5));  s = sin(ParG(5));
    a = ParG(3);  b = ParG(4);
    Xc = ParG(1);  Yc = ParG(2);
    P = (c/a)^2 - (s/b)^2;
    Q = (s/a)^2 - (c/b)^2;
    R = c*s*(1/a^2 + 1/b^2);
    ParA = [P; R; Q; -P*Xc-R*Yc; -Q*Yc-R*Xc; P*Xc^2+Q*Yc^2+2*R*Xc*Yc-1];
    ParA = ParA/norm(ParA);
end

if (code == 3)                           %   parabola
    c = cos(ParG(4));  s = sin(ParG(4));
    p = ParG(3);
    Xc = ParG(1);  Yc = ParG(2);
    R = Xc*s - Yc*c;
    ParA = [s^2; -c*s; c^2; -R*s-p*c; R*c-p*s; R^2+2*p*(Xc*c+Yc*s)];
    ParA = ParA/norm(ParA);
end

if (code == 4)                           %   intersecting lines
    c1=cos(ParG(1)); s1=sin(ParG(1));
    c2=cos(ParG(3)); s2=sin(ParG(3));
    ParA=[c1*c2; (c1*s2+c2*s1)/2;s1*s2; (c1*ParG(4)+c2*ParG(2))/2; (s1*ParG(4)+s2*ParG(2))/2; ParG(2)*ParG(4)];
    ParA = ParA/norm(ParA);
end

end
