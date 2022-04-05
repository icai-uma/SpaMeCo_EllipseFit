function [Handle]=PlotEllipseG(ParG,Color,linesize)
% Plot an ellipse, given its geometric parameters
if nargin < 3
    linesize = 1;
end
% Compute the natural parameters from the geometric parameters
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

% Obtain some points on the surface of the ellipse, in canonical
% coordinates, i.e. without rotation or shift
NumPoints=1000;
x1=zeros(1,NumPoints);
y1=zeros(1,NumPoints);
for i=1:NumPoints
    x1(i)=a*cos(2*pi*i/1000);
    y1(i)=b*sin(2*pi*i/1000);
end

% Rotate the points
V=[x1' y1']*[cos(phi) -sin(phi);sin(phi) cos(phi)];

% Shift the points
for i=1:1000
    x1(i)=V(i,1)+Center(1);
    y1(i)=V(i,2)+Center(2);
end

% Plot the points
Handle(1)=plot(x1,y1,'b-','Color',Color,'LineWidth',linesize);

% Plot the foci
Handle(2)=plot([Focus1(1) Focus2(1)],[Focus1(2) Focus2(2)],'bo','Color',Color,'LineWidth',linesize);
