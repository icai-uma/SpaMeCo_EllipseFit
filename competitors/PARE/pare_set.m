function x = pare_set (x, phi, alpha, a, b, z);
%PARE_SET       Param to vector conversion
%
% x = pare_set (x, phi, alpha, a, b, z);
% stores single parameters into param vector
%
% x: x == previous values of [phi; alpha; a; b; z]
% phi, alpha, a, b, z:
%   if (<val> != []), then set this param
%
% x: updated param vector

  m = size(x, 1) - 5;
  
  if (~isempty(phi)),   x(1:m) = phi; end;
  if (~isempty(alpha)), x(m+1) = alpha; end;
  if (~isempty(a)),     x(m+2) = a; end;
  if (~isempty(b)),     x(m+3) = b; end;
  if (~isempty(z)),     x(m+4:m+5) = z; end;

end % pare_set
