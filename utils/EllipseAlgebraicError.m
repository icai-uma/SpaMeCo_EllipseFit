function [Error]=EllipseAlgebraicError(TrueParA,ParA)
% Compute the error in the estimation of an ellipse, given the natural
% parametrization

% Parameter normalization
TrueParA = TrueParA/norm(TrueParA);
ParA = ParA/norm(ParA);

% Error computation
Error=norm(TrueParA-ParA);
% Error2=norm(TrueParA([3 4 1 2 5])-ParA);
% Error=min([Error1 Error2]);
