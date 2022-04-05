function [ParA,ParG,ParN,FinalError]=EllipseFitMunoz(X)
% Fit an ellipse to a set of training points.
% From:
% Muñoz-Pérez, J., et al. Multicriteria Robust Fitting of
% Elliptical Primitives. J Math Imaging Vis (2014) 49:492–509
% DOI 10.1007/s10851-013-0480-1
% Input:
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


% Algebraic method to initialize
%  A = EllipseFitByTaubin(X'); %Taubin (1991)
A = EllipseDirectFit(X'); %Fitzgibbon (1999)
% A = fit_ellipse(X(1,:), X(2,:)); %Halii y Flusser (1998)


u(1) = A(1);
u(2) = A(2);
u(3) = A(3);
u(4) = A(4);
u(5) = A(5);
u(6) = A(6);

show=1;
%u
[z, a, b, alpha, err] = ellipse_params(u, show);

%[z(1), z(2), a, b, alpha]
% %Geometric method (PARE algorithm)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meth=0;
% X=[x' y'];
% [z, a, b, alpha, phi, step] = pare (X, z, a, b, alpha, meth, show);
% 
% %%%%%%% End of the geometric method

% Computation of the foci
Center=z; % New definition of the ellipse center
c1=sqrt(abs(a^2-b^2));
Focus1=[Center(1)-c1*cos(alpha) Center(2)-c1*sin(alpha)];
Focus2=[Center(1)+c1*cos(alpha) Center(2)+c1*sin(alpha)];

  
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta=.01;       % Initial 
eta1=eta;       % step size for the gradient
eta2=eta;       %step size for the gradient
beta=0.01;      %penalty
factor1=1.05;    %1.05;  % Coefficient of increment of the step size
factor2=0.95;    %0.9;   % Coefficient of reduction of the step size
NumIterations=1000;    % Number of main loop iterations
error=zeros(1,NumIterations);



[~, NumPoints]=size(X);
x=X(1,:);
y=X(2,:);
d1=zeros(1,NumPoints);
d2=zeros(1,NumPoints);
d=zeros(1,NumPoints);
for i=1:NumPoints
     d1(i)=sqrt((x(i)-Focus1(1))^2+(y(i)-Focus1(2))^2);
     d2(i)=sqrt((x(i)-Focus2(1))^2+(y(i)-Focus2(2))^2);
     d(i)=d1(i)+d2(i);
end
m=median(d);

error(1)=sum(abs(d-m))/NumPoints + beta* sqrt((Focus1(1)-Focus2(1))^2 + (Focus1(2)-Focus2(2))^2); %NO METE EL BETA DE LA EQ. 13

% Variable initialization
xe1=0;
ye1=0;
xe2=0;
ye2=0;
xi1=0;
yi1=0;
xi2=0;
yi2=0;

nac=0;   
nev=0;
for k=2:NumIterations
    % Computation of the new foci
    for i=1:NumPoints
        if d(i)>m   % Exterior point
            xe1=xe1+(Focus1(1)-x(i))/(d1(i)+eps); % sum of the abscissae of the exterior points, x coordinate, 1st focus
            ye1=ye1+(Focus1(2)-y(i))/(d1(i)+eps); % sum of the ordinates of the exterior points, y coordinate, 1st focus
            xe2=xe2+(Focus2(1)-x(i))/(d2(i)+eps); % sum of the abscissae of the exterior points, x coordinate, 2nd focus
            ye2=ye2+(Focus2(2)-y(i))/(d2(i)+eps); % sum of the abscissae of the exterior points, y coordinate, 2nd focus
        else   % Interior point
            xi1=xi1+(Focus1(1)-x(i))/(d1(i)+eps); % sum of the abscissae of the interior points, x coordinate, 1st focus
            yi1=yi1+(Focus1(2)-y(i))/(d1(i)+eps); % sum of the ordinates of the interior points, y coordinate, 1st focus
            xi2=xi2+(Focus2(1)-x(i))/(d2(i)+eps); % sum of the abscissae of the interior points, x coordinate, 2nd focus
            yi2=yi2+(Focus2(2)-y(i))/(d2(i)+eps); % sum of the abscissae of the interior points, y coordinate, 2nd focus
        end
    end
    
    copiaf1=Focus1;
    copiaf2=Focus2;
    
    incrx1=eta1*(xe1-xi1);
    incry1=eta1*(ye1-yi1);
    incrx2=eta2*(xe2-xi2);
    incry2=eta2*(ye2-yi2);
    
    % Penalty term
    df = sqrt( (Focus1(1)-Focus2(1))^2 + (Focus1(2)-Focus2(2))^2 )+ eps;
    
    Focus1(1)=Focus1(1)-incrx1 - eta1*beta * ( Focus1(1)-Focus2(1) )/df;
    Focus1(2)=Focus1(2)-incry1 - eta2*beta * ( Focus1(2)-Focus2(2) )/df;
    Focus2(1)=Focus2(1)-incrx2 - eta1*beta * ( Focus2(1)-Focus1(1) )/df;
    Focus2(2)=Focus2(2)-incry2 - eta2*beta * ( Focus2(2)-Focus1(2) )/df;
    
   
    % Computation of the new distances and the error
    for i=1:NumPoints
        d1(i)=sqrt((x(i)-Focus1(1))^2+(y(i)-Focus1(2))^2);
        d2(i)=sqrt((x(i)-Focus2(1))^2+(y(i)-Focus2(2))^2);
        d(i)=d1(i)+d2(i);
    end
    m1=m;
    m=median(d);
    
    error(k)=sum(abs(d-m))/NumPoints + beta* sqrt((Focus1(1)-Focus2(1))^2 + (Focus1(2)-Focus2(2))^2); %NO METE EL BETA DE LA EQ. 13
    
    % Parameter update
    if error(k)<error(k-1)
        eta1=eta1*factor1;
        eta2=eta2*factor1;
        nac=nac+1; % Number of updates
    else
        eta1=eta1*factor2;
        eta2=eta2*factor2;
        Focus1=copiaf1;
        Focus2=copiaf2;
        m=m1;
        error(k)=error(k-1);
        nev=nev+1; % Number of evaluations
    end
    
end


% Compute the geometric parameters of the ellipse
Center=(Focus1+Focus2)/2;  % Ellipse center
a=m/2; % Length of the half major axis
b=0.5*sqrt(m^2-df^2); % Length of the half minor axis
phi=-atan(-(Focus2(2)-Focus1(2))/(Focus2(1)-Focus1(1))); % Tilt angle
ParG = [Center(1) Center(2) a b phi]';

% Compute the algebraic and natural parameters of the ellipse
ParA = GtoA(ParG,1);
ParN=[Focus1(1) Focus1(2) Focus2(1) Focus2(2) m]';

% Report the final error
FinalError=error(end);
