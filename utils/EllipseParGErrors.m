function Errors=EllipseParGErrors(TrueParG,ParG)
% Compute the errors in the estimation of an ellipse, given the geometric
% parametrization:
% Errors = [CenterError,AngleError,MajorSemiaxisError,MinorSemiaxisError,AreaError];


CenterError = norm(TrueParG(1:2)-ParG(1:2));

%  Making sure that ParG(3) is major, ParG(4) is minor (for true ellipse is sure)
if ParG(3) < ParG(4)  
    ab = ParG(3); ParG(3) = ParG(4); ParG(4) = ab;
    ParG(5) = ParG(5) - pi/2;
    if (ParG(5) < 0),  ParG(5) = ParG(5) + pi;  end
end

MajorSemiaxisError = abs(TrueParG(3)-ParG(3));
MinorSemiaxisError = abs(TrueParG(4)-ParG(4));

% making sure that bith angles are in first and second cuadrants.
if TrueParG(5) > pi, TrueParG(5) = TrueParG(5)-pi; end
if ParG(5) > pi, ParG(5) = ParG(5)-pi; end
AngleError = abs(TrueParG(5)-ParG(5));

AreaError = abs(pi*(TrueParG(3)*TrueParG(4)-ParG(3)*ParG(4)));

Errors = [CenterError,AngleError,MajorSemiaxisError,MinorSemiaxisError,AreaError];
