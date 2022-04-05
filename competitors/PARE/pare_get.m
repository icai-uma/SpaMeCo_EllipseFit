function [phi, alpha, a, b, z] = pare_get (x);
%PARE_GET       Vector to param conversion
%
% [phi, alpha, a, b, z] = pare_get (x);
% gets single parameters from param vector
%
% x: x == [phi; alpha; a; b; z]
%
% phi, alpha, a, b, z: splitted parameters

  m = size(x, 1) - 5;
  phi   = x(1:m); 
  alpha = x(m+1);
  a     = x(m+2);
  b     = x(m+3);
  z     = x(m+4:m+5);

end % pare_get
