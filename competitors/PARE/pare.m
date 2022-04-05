function [z, a, b, alpha, phi, step] = ...
         pare (X, z, a, b, alpha, meth, show);
%PARE   Geometric ellipse fit loop
%
% [z, a, b, alpha, phi, step] = ...
%      pare (X, z, a, b, alpha, meth, show{0});
% computes the best fit ellipse in parameterform
% x = z(1) + a  cos(phi-alpha), y = z(2) + b sin(phi-alpha)
%
% X: given points <X(i,1), X(i,2)>
% z, a, b, alpha: starting values for ellipse
% meth:
%   0 --> gauss-newton with marquardt for a near b
%   1 --> newton with marquardt for a near b
%   2 --> marquardt
%   3 --> gauss-newton
% show: if (show == 1), test output
%
% z, a, b, alpha: ellipse found
% phi: values for the nearest points (in parametric form)
% step: nof iterations
  
  if (nargin < 7), show = 0; end;
  if (nargin < 6), meth = 1; end;

  epsr  = 1e-5;
  phi=0; %Puesto por Oscar
  m = size (X, 1);
  x = zeros (m+5, 1);
  x = pare_set (x, phi, alpha, a, b, z);
  x = pare_initphi (X, x);

  step   = 0;
  normr  = 1;
  norma  = 1;
  lambda = [1;1;1];

  while (normr > epsr*norma),

    if     (meth == 0), [x, lambda] = pare_gauss_step (X, x, lambda);
    elseif (meth == 1), [x, lambda] = pare_newton_step (X, x, lambda);
    elseif (meth == 2), [x, lambda] = pare_marq_step (X, x, lambda);
    elseif (meth == 3), [x, lambda] = pare_gauss_step (X, x, [0;0;0]);
    else                error ('unknown meth');
    end

    if (step > 0),
      normr = norm (x - prevx);
      norma = norm (x);
    end

    prevx = x;
    step  = step+1;

    if (show == 1),
      [phi, alpha, a, b, z] = pare_get (x);
      drawellipse(z, a, b, alpha)
    end
    
    if (step > 100),
        if (show == 1),
            disp ('warning: number of steps exceeded limit');
        end
        break;
    end

  end % while

  [phi, alpha, a, b, z] = pare_get (x);
  
  %Añadido por Karl
  if a < b
    aux=a;
    a=b;
    b=aux;
    alpha=alpha+pi/2;
  end

end % pare
