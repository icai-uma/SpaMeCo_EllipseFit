function res = pare_residual (X, x);
%PARE_RESIDUAL  Residual vector for ellipse
%
% res = pare_residual (X, x);
% compute residual vector in the transformed system
% for parameter ellipse.
%
% X: given points <X(i,1), X(i,2)>
% x: ellipse parameters
%
% res: residual "X - ellipse" for the transformed system.
%   <res(i), res(m+i)> is residual for i-th point (m == nofpoints).
  
  [phi, alpha, a, b, z] = pare_get (x);

  s = sin(alpha);
  c = cos(alpha);
  Q = [c -s; s c];

  Xs = X*Q;
  zs = Q'*z;
  res = [Xs(:,1)-zs(1)-a*cos(phi); Xs(:,2)-zs(2)-b*sin(phi)];

end % pare_residual
