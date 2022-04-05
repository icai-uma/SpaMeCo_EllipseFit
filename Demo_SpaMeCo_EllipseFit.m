%% Demo of ellipse fitting with noise SpaMeCo_EllipseFit
%
% Please, cite this work as:
% 
% Karl Thurnhofer-Hemsi, Ezequiel López-Rubio, Elidia Beatriz Blázquez-Parra, M. Carmen Ladrón-de-Guevara-Muñoz, Óscar David de-Cózar-Macías,
% Ensemble ellipse fitting by spatial median consensus,
% Information Sciences, Volume 579, 2021, Pages 310-324, ISSN 0020-0255,
% https://doi.org/10.1016/j.ins.2021.08.011.
% (https://www.sciencedirect.com/science/article/pii/S0020025521008033)
%
% Last modification: 18/01/2021

clear all
warning off

% Prepare paths
addpath(genpath('./competitors'));
addpath(genpath('./utils'));
addpath('./L1mediancov/');

% Set seed
rng(8);

% Save results flag
saveV = false;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA GENERATION (SIMPLE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of ellipse generation
NumTestSamples = 1000;      % Number of test samples on the true ellipse used for testing
NumTrainSamples=50;         % Number of training samples
NoiseLevel=0.05;            % Noise level
OutlierProbability=0;    % Probability of outliers
OcclusionLevel=NaN;           % Level of occlusion
scale = 1;                  % Scale of the points (default in [0 1]x[0 1])
[TX,X,TrueParA,TrueParG,TrueParN]=GenerateRandomTestTrainingEllipse(NumTestSamples,NumTrainSamples,NoiseLevel,OutlierProbability,OcclusionLevel,scale);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN FITTING METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpaMeCo_EllipseFit
SubsamplingFactor=0.15; % Proportion of the total training samples considered for training
EnsembleSize=90;        % Number of subsamplings
percentile=10;          % Percentage(%) of best fits for the final ensemble
[ParA1,ParG1,ParN1]=SpaMeCo_EllipseFit(X);


%% Comparison with some methods
% Munoz (2014)
[ParA2,ParG2,ParN2]=EllipseFitMunoz(X);

% Fitzgibbon (1999)
A = EllipseDirectFit(X'); 
ParA3=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
[ParG3,code3]=AtoG(ParA3);
ParN3=GtoN(ParG3);

% Taubin (1991)
A = EllipseFitByTaubin(X');
ParA4=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
[ParG4,code4]=AtoG(ParA4);
ParN4=GtoN(ParG4);
if code4~=1
    ParG4 = NaN(5,1);
    ParA4 = NaN(6,1);
    ParN4 = NaN(5,1);
end

% Halii and Flusser (1998)
A = EllipseFitHalir(X(1,:), X(2,:)); 
ParA5=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
ParG5=AtoG(ParA5);
ParN5=GtoN(ParG5);

% Geometric method (PARE algorithm), initialization by Fitzgibbon (1999)
meth=0;
show=0;
[z, a, b, alpha, phi, step] = pare(X', ParG3(1:2), ParG3(3), ParG3(4), ParG3(5), meth, show);
ParG6=[z(1) z(2) a b alpha]';
ParA6=GtoA(ParG6,1);
ParN6=GtoN(ParG6);
if ~isreal(ParN6) || a<0 || b<0
    ParG6 = NaN(5,1);
    ParA6 = NaN(6,1);
    ParN6 = NaN(5,1);
end

% Rosin method, normalization A+C=1
[ParA7,ParG7,ParN7,code7]=EllipseFitRosin(X);
if (code7~=1)
    ParG7 = NaN(5,1);
    ParA7 = NaN(6,1);
    ParN7 = NaN(5,1);
end


% Szpak method
A = fast_guaranteed_ellipse_estimate(X');
ParA8 = [A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
[ParG8,code8]=AtoG(ParA8);
if code8==1
    ParN8=GtoN(ParG8); 
else
    ParG8 = NaN(5,1);
    ParA8 = NaN(6,1);
    ParN8=NaN(5,1);
end

% Prasad method
scaleUp = 100;
[residue,ea,eb,ex0,ey0,ealpha,code9] = Prasad_ElliFit(X(1,:)*scaleUp, X(2,:)*scaleUp);
if (code9==0)
    ParG9 = NaN(5,1);
    ParA9 = NaN(6,1);
    ParN9 = NaN(5,1);
else
    ParG9 = [ex0/scaleUp ey0/scaleUp ea/scaleUp eb/scaleUp ealpha]';
    ParA9=GtoA(ParG9,1);
    ParN9=GtoN(ParG9);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf,'Units','normalized');
set(gcf,'Position',[0.25,0.15,0.50,0.70]);

% Select a palette of colors
hold on
MyColors=distinguishable_colors(10);
MyColors(4,:)=MyColors(4,:)+0.5;
% Make yellow as True Ellipse
aux = MyColors(6,:);
MyColors(6,:) = MyColors(1,:);
MyColors(1,:) = aux;
% Switch color orders for the sake of clarity
aux = MyColors(6,:);
MyColors(6,:) = MyColors(5,:);
MyColors(5,:) = aux;

% Plot the ellipse, given the algebraic parameters
Handles=zeros(1,7);
if exist('TrueParG','var')
    MyHandle=PlotEllipseG(TrueParG,MyColors(1,:),4);
    Handles(1)=MyHandle(1);
end
MyHandle=PlotEllipseG(ParG2,MyColors(3,:));
Handles(3)=MyHandle(1);
MyHandle=PlotEllipseG(ParG3,MyColors(4,:));
Handles(4)=MyHandle(1);
MyHandle=PlotEllipseG(ParG4,MyColors(5,:));
Handles(5)=MyHandle(1);
MyHandle=PlotEllipseG(ParG5,MyColors(6,:));
Handles(6)=MyHandle(1);
MyHandle=PlotEllipseG(ParG6,MyColors(7,:));
Handles(7)=MyHandle(1);
MyHandle=PlotEllipseG(ParG7,MyColors(8,:));
Handles(8)=MyHandle(1);
MyHandle=PlotEllipseG(ParG8,MyColors(9,:));
Handles(9)=MyHandle(1);
MyHandle=PlotEllipseG(ParG9,MyColors(10,:));
Handles(10)=MyHandle(1);
% SAREfit is drawn at last
MyHandle=PlotEllipseG(ParG1,MyColors(2,:));
Handles(2)=MyHandle(1);

% Plot Training Samples
hold on
plot(X(1,:), X(2, :), '.r', 'Color', [0 0 0], 'MarkerSize',10);
xlabel('x');
ylabel('y');
axis equal
axis([-0.2 1.8 0 2])
grid on

% Save plot
if saveV
    PdfFileName=sprintf('./Demo_synthetic_rng%s',num2str(rngNum));
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperSize',[8 7]);
    set(gcf,'PaperPosition',[0 0 8 7]);
    set(gca,'fontsize',10);
    saveas(gcf,PdfFileName,'pdf');
else
    legend(Handles,'True','SpaMeCo_EllipseFit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','PARE','Rosin','Szpak','Prasad','Location','southoutside','Orientation','horizontal');
end

% Compute Natural errors
fprintf('Error for the new algorithm: %f\n',EllipseNaturalError(TrueParN,ParN1));
fprintf('Error for Munoz (2014) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN2));
fprintf('Error for Fitzgibbon (1999) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN3));
fprintf('Error for Taubin (1991) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN4));
fprintf('Error for Halii and Flusser (1998) algorithm: %f\n',EllipseNaturalError(TrueParN,ParN5));
fprintf('Error for PARE algorithm: %f\n',EllipseNaturalError(TrueParN,ParN6));
fprintf('Error for Rosin algorithm: %f\n',EllipseNaturalError(TrueParN,ParN7));
fprintf('Error for Szpak algorithm: %f\n',EllipseNaturalError(TrueParN,ParN8));
fprintf('Error for Prasad algorithm: %f\n',EllipseNaturalError(TrueParN,ParN9));


% Compute RMS Orthogonal Error
RSMO = zeros(1,9);

[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG1);
RSMO(1) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for the new algorithm: %f\n',RSMO(1));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG2);
RSMO(2) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Munoz (2014) algorithm: %f\n',RSMO(2));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG3);
RSMO(3) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Fitzgibbon (1999) algorithm: %f\n',RSMO(3));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG4);
RSMO(4) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Taubin (1991) algorithm: %f\n',RSMO(4));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG5);
RSMO(5) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Halii and Flusser (1998) algorithm: %f\n',RSMO(5));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG6);
RSMO(6) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for PARE algorithm: %f\n',RSMO(6));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG7);
RSMO(7) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Rosin algorithm: %f\n',RSMO(7));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG8);
RSMO(8) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Szpak algorithm: %f\n',RSMO(8));
[XYproj,~] = ProjectPointsOntoEllipse(TX',ParG9);
RSMO(9) = sqrt(mean(sqrt(sum((TX-XYproj).^2,1))));
fprintf('RSMO Error for Prasad algorithm: %f\n',RSMO(9));

