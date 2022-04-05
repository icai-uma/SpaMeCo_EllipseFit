function [z, a, b, alpha, err] = ellipse_params (u, show)
%ELLIPSE_PARAMS Get ellipse params from algebraic equation
%
%       [z, a, b, alpha, err] = ellipse_params (u, show{0});
%       get the ellipse parameters
%       from algebraic equation 
%         u(1)x^2 + u(2)xy + u(3)y^2 + ...
%         u(4)x + u(5)y + u(6) = 0.
%
%       u: coefficients of algebraic equation
%       show: == 1, then plot figure if error.
%
%       z, a, b, alpha: ellipse parameters
%       err: != 0, if not an ellipse

  if (nargin < 2) show = 0; end;
  err = 0;

  if u(1)<0
      u=-u;
  end
    
  A   = [u(1) u(2)/2; u(2)/2 u(3)];
  bb  = [u(4); u(5)]; 
  c   = u(6);

  [Q D] = eig(A);
  det   = D(1,1)*D(2,2);
  if (det <= 0),
    err = 1;
    if (show == 1), drawconic (u); end;
    z = [0;0];
    a = 1; b = 1; alpha = 0;
  else 
    bs    = Q'*bb;
    alpha = atan2(Q(2,1), Q(1,1));
    zs    = -(2*D)\bs;  
    z     = Q*zs;
    h     = -bs'*zs/2-c;
    a     = sqrt(h/D(1,1));
    b     = sqrt(h/D(2,2));
  end

end % ellipse_params  
