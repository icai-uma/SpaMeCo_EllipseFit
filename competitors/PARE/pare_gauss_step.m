function [x, lambda] = pare_gauss_step(X, x, lambda);
%PARE_GAUSS_STEP        Gauss-Newton iteration step
%
% [x, lambda] = pare_gauss_step(X, x, lambda);
% makes basic step for this x. Adds marquart correction if 
% abs(a - b)/(a + b) < lambda(1).
%
% X: given points <X(i,1), X(i,2)>
% x: given parameters
% lambda: lambda(1) marquardt correction factor
%
% x: updated parameters
% lambda: (possibly new) marquardt factor

  [phi, alpha, a, b, z] = pare_get (x);

  if (abs(a - b)/(a + b) < lambda(1)), 
    [x, lambda] = pare_marq_step(X, x, lambda);
  else
    m = size(X,1);
    Y = pare_residual (X, x);

  %% form Jacobian
    S = sin(phi);  
    C = cos(phi);
    s = sin(alpha);
    c = cos(alpha);
    A = [-b*S C zeros(size(phi))  c*ones(size(phi)) s*ones(size(phi))];
    B = [ a*C zeros(size(phi)) S -s*ones(size(phi)) c*ones(size(phi))];

    [cg, sg] = rot_cossin (-a*S, b*C);
    G = sparse ([diag(cg), -diag(sg); diag(sg), diag(cg)]);
    Y = G*Y;

    D = diag (- a*S.*cg - b*C.*sg);
    J = [[D; zeros(m,m)], G*[A; B]];

    [qq, J(m+1:2*m, m+1:m+5)] = qr(J(m+1:2*m, m+1:m+5));    
    Y(m+1:m+5, :) = qq(:,1:5)'*Y(m+1:2*m,:);

    h = J(1:m+5, 1:m+5)\Y(1:m+5);

    x = x + h;
     
  end % if (marq term necessary)

end % pare_gauss_step
