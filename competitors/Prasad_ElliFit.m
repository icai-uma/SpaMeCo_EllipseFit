% Usage: [residue,a,b,x0,y0,alpha,status]=Prasad_ElliFit(xt,yt);
%
% Inputs:  
% xt - the x-coordinates of the data points
% yt - the y-coordinates of the data points
% minimum no of (xt,yt) pairs should be greater than equal to 7
%    
% Outputs:
% residue: least squares residue
% a: semi-major axis
% b: semi-minor axis
% x0: x-ccordinate of the center of the ellipse
% y0: y-coordinate of the center of the ellipse
% alpha: the orientation angle of the ellipse (please check - based on your
%        definition of the coordinates, you may need to do alpha=-alpha)
% status: 0 if ellipse cannot be fit, 1 otherwise
%     
% This function fits ellipse on a set of pixels (supposedly on the
% edge/boundary) using geometry based unconstrained, non-iterative least
% squares method.
%
% Copyright (c) 2012 Dilip K. Prasad
% School of Computer Engineering
% Nanyang Technological University, Singapore
% http://www.ntu.edu.sg/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software with restriction for its use for research purpose only,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% Please cite the following work if this program is used:
% [1] Dilip K. Prasad, Maylor K.H. Leung and Chai Quek, "ElliFit: An 
% unconstrained, non-iterative, least squares based geometric Ellipse 
% Fitting method,”  Pattern Recognition, vol. 46 issue 5,pp. 1449-1465, 2013.

function [residue,ea,eb,ex0,ey0,ealpha,estatus]=Prasad_ElliFit(xt,yt)
warning off;
xt_act=xt(:);
yt_act=yt(:);
xt_shift=(min(xt_act)+max(xt_act))/2;
yt_shift=(min(yt_act)+max(yt_act))/2;
xt=xt(:)-xt_shift;yt=yt(:)-yt_shift;

A=[(xt).^2 2*xt.*yt -2*xt -2.*yt -1*ones(size(xt))];
B=-yt.^2;
Phi=inv(A.'*A)*A.'*B;

residue=norm(A*Phi-B)/norm(B);
estatus=1;
[ea,eb,ex0,ey0,ealpha]=compute_ellipse_parameter_from_Phi(Phi);
ex0=ex0+xt_shift;
ey0=ey0+yt_shift;
ealpha=-ealpha;

if ~isreal(ea) || ~isreal(eb) || isnan(ea) || isnan(eb)
    estatus=0;
    ea=0;eb=0;ex0=0;ey0=0;ealpha=0;
end
end

function [a,b,x0,y0,alpha]=compute_ellipse_parameter_from_Phi(Phi)
x0=(Phi(3)-Phi(4)*Phi(2))./((Phi(1))-(Phi(2))^2);
y0=(Phi(1)*Phi(4)-Phi(3)*Phi(2))./((Phi(1))-(Phi(2))^2);
term2=sqrt(((1-Phi(1))^2+4*(Phi(2))^2));
term3=(Phi(5)+(y0)^2+(x0^2)*Phi(1)+2*Phi(2));
term1=1+Phi(1);
b=(sqrt(2*term3/(term1+term2)));
a=(sqrt(2*term3/(term1-term2)));
alpha=0.5*atan2(2*Phi(2),1-Phi(1));

% x0=(Phi(2)*Phi(4)-2*Phi(3))./((Phi(2))^2-4*Phi(1));
% y0=(Phi(2)*Phi(3)-2*Phi(1)*Phi(4))./((Phi(2))^2-4*Phi(1));
% term2=sqrt(((1-Phi(1))^2+(Phi(2))^2));
% term3=(Phi(2)*Phi(3)*Phi(4)-Phi(4)^2*Phi(1)-Phi(3)^2)/(Phi(2)^2-4*Phi(1))-Phi(5);
% term1=1+Phi(1);
% b=2*(sqrt(term3/(term1+term2)));
% a=2*(sqrt(term3/(term1-term2)));
% alpha=-0.5*atan2(-2*Phi(2),1-Phi(1));

end
