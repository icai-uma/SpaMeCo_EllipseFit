% Demo of the Euclidean Ellipse Comparison Metric (EECM),
% an ellipse fit error and ellipse distance calculation method:
%  
% Cakir, H. I., Topal, C. "An Euclidean Ellipse Comparison Metric 
% for Quantitative Evaluation", International Conference on 
% Acoustics, Speech and Signal Processing (ICASSP), 2018.

clear all
close all

MyEllipse=[1 2 3 4 5];
% MyEllipse2=[11 22 7 40 5];
MyEllipse2=[1 2 3 4 5];
MyPixel=[3.4 2.5];

noPoints=100;
points=CalculateEllipsePointsMEX(MyEllipse,noPoints);
points2=CalculateEllipsePointsMEX(MyEllipse2,noPoints);
hold on
plot(points(1,:),points(2,:),'*g')
plot(points2(1,:),points2(2,:),'*r')
% plot(MyPixel(1),MyPixel(2),'*b')
plot(points(1,3),points(2,3),'*b')
axis equal
hold off

MyDistance=DistanceToPointMEX(MyEllipse, MyPixel)

MyDistance0=DistanceToPointMEX(MyEllipse, points(:,3))

MyDistanceEllipses=DistanceBwEllipsesMEX(MyEllipse, MyEllipse2, noPoints)
