function x = pare_initphi (X, x);
%PARE_INITPHI
%
%       x = pare_initphi (X, x);
%       assigns approximate values for angles
%       relative to the transformed system.
%
%       X: given points <X(i,1), X(i,2)>
%       x: parameters
%
%       x: parameters, with 'phi' values approximate
%       nearest point angles

  [phi, alpha, a, b, z] = pare_get (x);

  c = cos(alpha); s = sin(alpha);
  Q = [c -s; s c];
  % compute initial approximations for phi_i
  du  = Q'*[ X(:,1)-z(1) X(:,2)-z(2)]';
  phi = angle(du(1,:)/a  + sqrt(-1)*du(2,:)/b)';

  x = pare_set (x, phi, [], [], [], []);

end % pare_initphi
