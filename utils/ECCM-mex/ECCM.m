function [Error]=ECCM(ParG1,ParG2,nPoints)
% Compute the Euclidean Ellipse Comparison Metric (EECM),
% an ellipse fit error and ellipse distance calculation method:
%  
% Cakir, H. I., Topal, C. "An Euclidean Ellipse Comparison Metric 
% for Quantitative Evaluation", International Conference on 
% Acoustics, Speech and Signal Processing (ICASSP), 2018.

if any(isnan(ParG1)) || any(isnan(ParG2))
    Error = NaN;
else
    ParG1(5) = rad2deg(ParG1(5));
    ParG2(5) = rad2deg(ParG2(5));
    Error = DistanceBwEllipsesMEX(ParG1,ParG2,nPoints);
end



end