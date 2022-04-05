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
rng('default');

% Path of the images
ImagePath = './images';
%% DataEllipse
Images = {'Hda_obj93','7bg','saturn'};
Suffix = {'','','',''};


%% Batch running all images
NumImages=numel(Images);
NumMethods = 9;
Labels = {'SpaMeCo_EllipseFit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','PARE','Rosin','Szpak','Prasad'};

for NdxImage=1:NumImages
    
    ThisImage = Images{NdxImage};
    ThisSuffix = Suffix{NdxImage};
    ThisData = [ThisImage ThisSuffix];

    % Load MAT file with detected points (done with the desired algorithm)
    % It should have the same name of the image and suffix 'p'
    MatImage = strcat(sprintf('%s/p',ImagePath),ThisData);
    load(MatImage)
    X = [y;x];

    % Execute the methods
        
    % SpaMeCo_EllipseFit
    [ParA1,ParG1,ParN1]=SpaMeCo_EllipseFit(X);

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
    [z, a, b, alphaA, phi, step] = pare(X', ParG3(1:2), ParG3(3), ParG3(4), ParG3(5), meth, show);
    ParG6=[z(1) z(2) a b alphaA]';
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

    
    % Plot results
    MyColors=distinguishable_colors(NumMethods+1);
    
    % Make grey the black one
    MyColors(4,:)=MyColors(4,:)+0.5;
    
    % Make yellow as True Ellipse
    aux = MyColors(6,:);
    MyColors(6,:) = MyColors(1,:);
    MyColors(1,:) = aux;
    
    % Switch color orders for the sake of clarity
    aux = MyColors(6,:);
    MyColors(6,:) = MyColors(5,:);
    MyColors(5,:) = aux;

    % Plot image with points and ellipses
    figure
    JpgImage = imread(strcat(sprintf('%s/',ImagePath),ThisImage,'.jpg'));
    imshow(JpgImage)
    PositionO = get(gca,'OuterPosition');
    PositionI = get(gca,'InnerPosition');
    alpha 0.7
    hold on
    Handles=zeros(1,NumMethods+1);
    if exist('TrueParG','var')
        MyHandle=PlotEllipseG(TrueParG,MyColors(1,:),2);
        Handles(1)=MyHandle(1);
    end
    MyHandle=PlotEllipseG(ParG2,MyColors(3,:),0.75);
    Handles(3)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG3,MyColors(4,:),0.75);
    Handles(4)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG4,MyColors(5,:),0.75);
    Handles(5)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG5,MyColors(6,:),0.75);
    Handles(6)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG6,MyColors(7,:),0.75);
    Handles(7)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG7,MyColors(8,:),0.75);
    Handles(8)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG8,MyColors(9,:),0.75);
    Handles(9)=MyHandle(1);
    MyHandle=PlotEllipseG(ParG9,MyColors(10,:),0.75);
    Handles(10)=MyHandle(1);
    % SAREfit is drawn at last
    MyHandle=PlotEllipseG(ParG1,MyColors(2,:),0.75);
    Handles(2)=MyHandle(1);

    % Draw samples
    plot(X(1,:), X(2, :), '+r', 'Color', MyColors(1,:),'MarkerSize',2);
    axis equal tight
    axis([0 size(JpgImage,2) 0 size(JpgImage,1)])
    set(gca,'YDir','reverse')
    set(gca,'OuterPosition',PositionO)
    set(gca,'InnerPosition',PositionI)

    PdfFileName=sprintf('%s/%s_ImageEllipses',ImagePath,ThisData);
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperSize',[8 7]);
    set(gcf,'PaperPosition',[0 0 8 7]);
    set(gca,'fontsize',10);
    saveas(gcf,PdfFileName,'pdf');
   
%     legend(Handles,'SpaMeCo_EllipseFit','Muñoz','Fitzgibbon','Taubin','Halir&Flusser','PARE','Rosin','Szpak','Prasad','Location','southoutside','Orientation','horizontal');
    
end




