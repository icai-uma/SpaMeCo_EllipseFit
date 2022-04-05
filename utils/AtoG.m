function [ParG, code] = AtoG (ParA)
% http://people.cas.uab.edu/~mosya/cl/MATLABconics.html

%  Conversion of Algebraic parameters of a conic to its Geometric parameters

%   Algebraic parameters are coefficients A,B,C,D,E,F in the algebraic
%   equation     Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0

%   Geometric parameters depend on the type of the conic (ellipse, etc.)

%   Input:  ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters

%   Output: code is the code of the conic type (see below)
%           ParG is the vector of Geometric parameters (see below)

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

ParA=ParA/norm(ParA);    %  normalize the given algebraic parameter vector

ParG = -1;   %  this will be returned for imaginary or degenerate conics

if (1+ParA(1)==1) && (1+ParA(2)==1) && (1+ParA(3)==1)  %  if no quadratic part...
    if (1+ParA(4)==1) && (1+ParA(5)==1)
        code = 11;     %  the "pole", extreme singularity
    else
        code = 7;      %  single line (no quadratic part, only linear part)
    end
    return;
end

M33 = [ParA(1) ParA(2) ParA(4);
       ParA(2) ParA(3) ParA(5);
       ParA(4) ParA(5) ParA(6)];   %    big, 3x3 matrix
M22 = [ParA(1) ParA(2);
       ParA(2) ParA(3)];           %  small, 2x2 matrix

det3x3 = det(M33);                 %    3x3 determinant
det2x2 = det(M22);                 %    2x2 determinant

if (1+det3x3 == 1)                   %  if the big matrix is singular...
    if (1+det2x2 == 1)
        dettwo = ParA(1)*ParA(6) - ParA(4)^2 + ParA(3)*ParA(6) - ParA(5)^2;
        if (dettwo > 0),  code = 10;  end   %  imaginary parallel lines
        if (dettwo < 0),  code =  5;  end   %  parallel lines
        if (1+dettwo == 1),  code = 6;  end   %  coincident lines
        return;
    end
    if (det2x2 > 0)
        code = 8;                %  single point
    else
        code = 4;                            %  intersecting lines
        [Q D] = eig(M33);    %  eigendecomposition
        [Dmax,Imax] = max(diag(D));
        [Dmin,Imin] = min(diag(D));
        Qmax = Q(1:3,Imax)*sqrt(abs(Dmax));
        Qmin = Q(1:3,Imin)*sqrt(abs(Dmin));
        Q1 = Qmax + Qmin;
        Q2 = Qmax - Qmin;
        theta1 = atan2(Q1(2),Q1(1));
        d1 = Q1(3)/norm(Q1(1:2));
        theta2 = atan2(Q2(2),Q2(1));
        d2 = Q2(3)/norm(Q2(1:2));
        ParG = [theta1; d1; theta2; d2];
    end
    return;
end

%          Next: non-degenrate types of conics

[Q D] = eig(M22);    %  eigendecomposition

U  = Q'*[ParA(4); ParA(5)];       %  orthogonal transformation

if (1+D(1,1)==1)||(1+D(2,2)==1)
    code = 3;                             %   parabola
    if (abs(D(1,1))>abs(D(2,2)))
        Uc1 = -U(1)/D(1,1);
        Uc2 = -(U(1)*Uc1 + ParA(6))/2/U(2);
        Center = Q*[Uc1; Uc2];
        p = -U(2)/D(1,1);
        Angle = atan2(Q(2,2), Q(1,2));
        ParG = [Center; p; Angle];
    else
        Uc2 = -U(2)/D(2,2);
        Uc1 = -(U(2)*Uc2 + ParA(6))/2/U(1);
        Center = Q*[Uc1; Uc2];
        p = -U(1)/D(2,2);
        Angle = atan2(Q(2,1), Q(1,1));
        ParG = [Center; p; Angle];
    end
    return;
end

Uc = -U./diag(D);
Center = Q*Uc;
H = -U'*Uc - ParA(6);
if (D(1,1)*D(2,2) < 0)
    code = 2;                            %   hyperbola
    if (D(1,1)*H > 0)
        a = sqrt(H/D(1,1));
        b = sqrt(-H/D(2,2));
        Angle = atan2(Q(2,1), Q(1,1));
        if (Angle < 0),  Angle = Angle + pi;    end
    else
        a = sqrt(H/D(2,2));
        b = sqrt(-H/D(1,1));
        Angle = atan2(Q(2,2), Q(1,2));
        if (Angle < 0),  Angle = Angle + pi;    end
    end
    ParG = [Center; a; b; Angle];
else
    if (H*D(1,1) <= 0)
        code = 9;                     %  imaginary ellipse
    else
        code = 1;                        %   ellipse
        a = sqrt(H/D(1,1));
        b = sqrt(H/D(2,2));
        Angle = atan2(Q(2,1), Q(1,1));
        if (Angle < 0),  Angle = Angle + pi;    end
        if a < b      %  making sure that a is major, b is minor
            ab = a; a = b; b = ab;
            Angle = Angle - pi/2;
            if (Angle < 0),  Angle = Angle + pi;  end
        end
        ParG = [Center; a; b; Angle];
    end
end
    
end    %  end of function AtoG
