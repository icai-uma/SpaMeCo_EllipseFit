function perf(T,metric,labels,logplot)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004

% load('./Datos/MSE_QR.mat')
% load('./Datos/MSE_wQR.mat')
% load('./Datos/MSE_R1wQR.mat')
% load('./Datos/MSE_iwQR2.mat')
% load('./Datos/MSE_SVD.mat')
% load('./Datos/MSE_HH.mat')
% load('./Datos/MSE_HHm.mat')
% load('./Datos/MSE_T.mat')
% load('./Datos/MSE_Ainv.mat')
% load('./Datos/MSE_HH_T.mat')
% load('./Datos/MSE_CSNE.mat')
% load('./Datos/MSE_CSNE_T.mat')
% load('.\Datos\CPUtime_QR.mat')
%     nivel = 1;
% T = [MSE_QR MSE_wQR MSE_R1wQR MSE_iwQR2...
%     MSE_SVD MSE_HH MSE_HHm MSE_HH_T...
%     MSE_T MSE_Ainv MSE_CSNE MSE_CSNE_T];


if (nargin < 3)
    logplot = 0; 
end
% logplot = 0;

% colors  = ['m' 'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b' 'r' 'g' 'k'];
colors = distinguishable_colors(10);

% Make grey the black one
colors(4,:)=colors(4,:)+0.5;

% Make yellow as True Ellipse
aux = colors(6,:);
colors(6,:) = colors(1,:);
colors(1,:) = aux;

% Switch color orders for the sake of clarity
aux = colors(6,:);
colors(6,:) = colors(5,:);
colors(5,:) = aux;

lines   = [':' '-' '-.' '--' ':' '-' '-.' '--' ':' '-' '-.' '--'];
markers = ['x' '*' 's' 'd' 'v' '^' 'o' 'v' '^' 'o' 'x' '*' 's' 'd'];

[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;
r = sort(r);

% Plot stair graphs with markers.

figure
for s = 1: ns
 [xs,ys] = stairs(r(:,s),[1:np]/np);
%  option = ['-' colors(s) markers(s)];
    plot(xs,ys,'Color',colors(s+1,:),'MarkerSize',14-s);
    hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

axis([ 0 1.1*max_ratio 0 1 ]);

% Legends and title should be added.
% legend(labels);
% title('Performance profile con el error como métrica y un problema por nivel');
xlabel(sprintf('\\tau    (%s)',metric));
if logplot, xlabel(sprintf('log_2(\\tau)    (%s)',metric)); end
ylabel('P(r_{p,m}\leq\tau : 1\leq m\leq n_m)');
