function [x, lambda] = pare_marq_step(X, x, lambda);
%PARE_MARQ_STEP         Gauss-Newton step with Marquardt
%
% [x, lambda] = pare_marq_step(X, x, lambda);
% makes basic step for this x with marquardt correction. 
%
% X: given points <X(i,1), X(i,2)>
% x: given parameters
% lambda: lambda(i) i-th previous marquardt correction factor
%
% x: updated parameters
% lambda: updated marquardt factors

  [phi, alpha, a, b, z] = pare_get (x);
  
  mu      = 1e-3;
  omega   = 0.5;

  W = ones(size(x,1), 1);
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
  G = sparse([diag(cg), -diag(sg); diag(sg), diag(cg)]);
  Y = G*Y;

  D  = - a*S.*cg - b*C.*sg;
  AB = G*[A; B];

%% do marquart step
  istep = 0;
  while (1),
    [cg, sg] = rot_cossin (D, lambda(1)*ones(m,1));
    DD = D.*cg - lambda(1)*sg;
    AA = sparse(diag(cg))*AB(1:m,:);
    YY = [sparse(diag(cg))*Y(1:m,:); Y(m+1:2*m,:); ...
          sparse(diag(sg))*Y(1:m,:)]; 
    BB = [AB(m+1:2*m,:); sparse(diag(sg))*AB(1:m,:); ...
          lambda(1)*eye(5,5)];
    [qq, RR] = qr(BB);
    YY = [YY(1:m,:); qq(1:2*m,1:5)'*YY(m+1:3*m,:)];
    JJ = [diag(DD), AA; zeros(5,m), RR(1:5,:)];
    s  = JJ\YY;

    h = norm(Y) - norm(pare_residual (X, x + s));
    if (h >= 0), break; end;
    lambda(1) = lambda(1)/omega;
    istep = istep + 1;
  end % while
  if (istep == 0),
    lambda(1) = lambda(1)*omega;
  end

  x = x + s;

end % pare_marq_step
