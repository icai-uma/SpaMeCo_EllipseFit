function [Error]=EllipseNaturalError(TrueParN,ParN)
% Compute the error in the estimation of an ellipse, given the natural
% parametrization
Error1=norm(TrueParN-ParN);
Error2=norm(TrueParN([3 4 1 2 5])-ParN);
Error=min([Error1 Error2]);
