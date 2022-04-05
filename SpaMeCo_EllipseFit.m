function [ParA,ParG,ParN]=SpaMeCo_EllipseFit(X)
% Fit an ellipse to a set of training points.
% Consensus of basic fits.
% Inputs:
%   X=matrix of size 2 x NumPoints with the training samples
% Outputs:
%   ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]' is the vector of 
%   geometric parameters of the ellipse). a=half length of the major axis,
%   b=half length of the minor axis
%   ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters:
%           Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0
%   ParN= [Focus1x Focus1y Focus2x Focus2y SumDists] is the vector of
%   natural parameters of the ellipse: the two foci (Focus1 and Focus2),
%   and the sum of distances to both foci, SumDists==2a


EnsembleSize=5;

% Vector to control methods which do not fit properly
error = 0;
NdxNumMethodsUsed = 1;

% Compute the members of the ensemble
for NdxEnsemble=1:EnsembleSize
    switch NdxEnsemble
        case 1
            % Taubin (1991)
            A = EllipseFitByTaubin(X');
            ThisParA=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
            [ThisParG,code]=AtoG(ThisParA);
            ThisParG = real(ThisParG);
            ThisParN=GtoN(ThisParG);
            if code~=1, error = 1; end

            
        case 2
            % Geometric method (PARE algorithm), initialization by Fitzgibbon (1999)
            A = EllipseDirectFit(X');
            ParAIni=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
            ParGIni=real(AtoG(ParAIni));
            meth=0;
            show=0;
            [z, a, b, alpha, ~, ~] = pare(X', ParGIni(1:2), ParGIni(3), ParGIni(4), ParGIni(5), meth, show);
            ThisParG=[z(1) z(2) a b alpha]';
            ThisParN=GtoN(ThisParG);
            if ~isreal(ThisParN) || a<0 || b<0, error = 1; end
            
        case 3
            % Fitzgibbon (1999)
            A = EllipseDirectFit(X');
            ThisParA=[A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
            ThisParG=real(AtoG(ThisParA));
            ThisParN=GtoN(ThisParG);

        case 4
            % Munoz method
            [~,ThisParG,ThisParN]=EllipseFitMunoz(X);
            LastResource = ThisParG;


        case 5
            % Szpak method
            A = fast_guaranteed_ellipse_estimate(X');
            ThisParA = [A(1) 0.5*A(2) A(3) 0.5*A(4) 0.5*A(5) A(6)]';
            [ThisParG,code]=AtoG(ThisParA);
            if code~=1, error=1; 
            else
                ThisParG = real(ThisParG);
                ThisParN=GtoN(ThisParG);
                LastResource = ThisParG; 
            end

    end

    
    if ~error
        % Store the geometric parameters of this ensemble member
        AllParG(:,NdxNumMethodsUsed)=ThisParG;

        % Obtain the natural parameters of the ellipse
        AllParN(:,NdxNumMethodsUsed)=ThisParN;
        
        NdxNumMethodsUsed = NdxNumMethodsUsed+1;
    end
    error = 0;

end

fprintf('Number of methods in the consensus: %d\n',NdxNumMethodsUsed-1)

% Compute the consensus from the natural parameters of the ensemble members
ParN=zeros(5,1);
ParN(1:2)=L1mediancovMEX(AllParN(1:2,:));
ParN(3:4)=L1mediancovMEX(AllParN(3:4,:));
ParN(5)=median(AllParN(5,:));
% Compute geometric parameters
Focus1=ParN(1:2);
Focus2=ParN(3:4);
Center=(Focus1+Focus2)/2;  % Ellipse center
df=norm(Focus1-Focus2); % Focal distance
m=ParN(5);

% Check whether this is an ellipse
if m>df
    % It is an ellipse, so it is OK
    a=m/2; % Length of the half major axis
    b=0.5*sqrt(m^2-df^2); % Length of the half minor axis
    phi=atan((Focus2(2)-Focus1(2))/(Focus2(1)-Focus1(1))); % Tilt angle
    ParG=[Center(1) Center(2) a b phi]';
else
    % It is not an ellipse, so we revert to the Muñoz solution
    ParG=LastResource;
    disp('ERROR: The consensus is not an ellipse. Using emergency solution')
end

% Convert to algebraic form
ParA=GtoA(ParG,1);

% Convert to natural form
ParN=GtoN(ParG);

