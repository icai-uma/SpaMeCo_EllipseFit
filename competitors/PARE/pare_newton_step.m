function [x, lambda] = pare_newton_step (X, x, lambda);
%PARE_NEWTON_STEP       Newton iteration step
%
% [x, lambda] = pare_newton_step (X, x, lambda);
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
    ALPHA = m + 1;
    A     = m + 2;
    B     = m + 3;
  
    Y = pare_residual (X, x);
    YY = [Y(1:m, 1), Y(m+1:2*m, 1)];
        
  %% form Jacobian
    S = sin(phi);  
    C = cos(phi);
    s = sin(alpha);
    c = cos(alpha);
    JA = [-b*S C zeros(size(phi))  c*ones(size(phi)) s*ones(size(phi))];
    JB = [ a*C zeros(size(phi)) S -s*ones(size(phi)) c*ones(size(phi))];

    DA = -a*sparse(diag(S));
    DB = b*sparse(diag(C));

    H  = zeros(m+5, m+5);
    for i = 1:m, H(i,i)     = [a*C(i), b*S(i)]*YY(i,:)'; end;
    for i = 1:m, H(i,ALPHA) = [b*C(i), a*S(i)]*YY(i,:)'; end;
    H(1:m,A)                = S.*YY(:,1);
    H(1:m,B)                = -C.*YY(:,2);
    H(ALPHA, ALPHA)         = [a*C', b*S']*Y; 
    H(ALPHA, A)             = C'*YY(:,2); 
    H(ALPHA, B)             = -S'*YY(:,1);
    H = H + triu(H, 1)';

    DD = a^2*S.^2 + b^2*C.^2;
    J1 = DA*JA + DB*JB;
    J2 = [diag(DD), J1; J1', JA'*JA + JB'*JB];
    Y2 = [DA*Y(1:m) + DB*Y(m+1:2*m); JA'*Y(1:m) + JB'*Y(m+1:2*m)];
    s = (J2 + H)\Y2;

    x = x + s;

  end % if (marq term necessary)

end % pare_newton_step
