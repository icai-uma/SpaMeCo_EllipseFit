function [c, s] = rot_cossin (x, y);
%ROT_COSSIN     Givens rotation angles
%
% [c, s] = rot_cossin (x, y);
% returns cos and sin vectors for Givens-rotation matrix
% which rotates y to zero.
%
% x, y: vectors
% c(i), s(i): [c(i) -s(i); s(i) c(i)]*[x(i); y(i)] == [..; 0]
  
  m = size(x,1);
  c = zeros(m,1); s = zeros(m,1);
  for i=1:m,
    if (abs(y(i)) > abs(x(i))),
      cot = -x(i)/y(i); si = 1/sqrt(1+cot^2); co = si*cot;
    else
      tan = -y(i)/x(i); co = 1/sqrt(1+tan^2); si = co*tan;
    end
    s(i) = si; c(i) = co;
  end 

end % rot_cossin
