function [TestSamples,TrainingSamples,ParA,ParG,ParN]=GenerateRandomTestTrainingEllipse(NumTestSamples,NumTrainSamples,NoiseLevel,OutlierProbability,OcclusionLevel,scale)
% Generate a random ellipse and draw some noisy samples from it
% Inputs:
%   NumSamples=The number of random samples to be drawn
%   NoiseLevel=positive real number which specifies the standard deviation
%   of the Gaussian noise to be added to the samples
% Outputs:
%   Samples=The generated noisy samples
%   ParG = [Xcenter, Ycenter, a, b, AngleOfTilt]' is the vector of 
%   geometric parameters of the ellipse). a=half length of the major axis,
%   b=half length of the minor axis
%   ParA = (A,B,C,D,E,F)' is the vector of Algebraic parameters:
%           Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0
%   ParN= [Focus1x Focus1y Focus2x Focus2y SumDists] is the vector of
%   natural parameters of the ellipse: the two foci (Focus1 and Focus2),
%   and the sum of distances to both foci, SumDists==2a


% Choose the center in the unit square [0,1]x[0,1]
Center=scale*rand(2,1);

% Choose the major and minor axes
a=scale*(0.2+0.8*rand(1));
b=scale*(0.1+0.9*rand(1));
% Ensure that the major axis is bigger than the minor axis
if a<b
    Temporary=a;
    a=b;
    b=Temporary;
end

% Choose the tilt angle
phi=pi*rand(1)-pi/2;

% Compute some parameters of the ellipse
% CenterToFocusDistance=sqrt(a^2-b^2);
% Focus1=[-CenterToFocusDistance 0]';
% Focus2=[CenterToFocusDistance 0]';
% RotationMatrix=[cos(-phi) -sin(-phi);sin(-phi) cos(-phi)];
% Focus1=RotationMatrix*Focus1+Center;
% Focus2=RotationMatrix*Focus2+Center;



% Choose a starting parameter and an ending parameter on the ellipse to generate
% the samples on that interval. We want an interval which is larger than
% 1, in order to avoid datasets with too small curvature which lead to
% degenerate solutions.
Done=0;
while ~Done
    StartingParameter=2*pi*rand(1)-pi;
    if isnan(OcclusionLevel)    % Random occlusion
        EndingParameter=2*pi*rand(1)-pi;
    else
        EndingParameter=StartingParameter+2*pi*(1-OcclusionLevel);
    end
    if (StartingParameter>EndingParameter)
        Temporary=StartingParameter;
        StartingParameter=EndingParameter;
        EndingParameter=Temporary;
    end
    Done=(EndingParameter-StartingParameter)>1;     % To choose the lenght of the arc (in radians) ->  Antes era 1 (o pi/3)
end

% Generate the samples on the canonical coordinate system, i.e. without
% rotation or shift
Samples=zeros(2,NumTrainSamples);
for NdxSample=1:NumTrainSamples
    MyAngle=StartingParameter+(EndingParameter-StartingParameter)*rand(1);
    Samples(1,NdxSample)=a*cos(MyAngle);
    Samples(2,NdxSample)=b*sin(MyAngle); 
end

% Add Gaussian noise
TrainingSamples=Samples+NoiseLevel*scale*randn(2,NumTrainSamples);

% Add outliers in the box [-1,1] x [-1,1]
OutlierIndices=find(rand(1,NumTrainSamples)<OutlierProbability);
TrainingSamples(:,OutlierIndices)=-1+2*rand(2,numel(OutlierIndices));

% Generate test samples dataset
TestSamples=zeros(2,NumTestSamples);
for NdxSample=1:NumTestSamples
    MyAngle=-pi+2*pi*rand(1);
    TestSamples(1,NdxSample)=a*cos(MyAngle);
    TestSamples(2,NdxSample)=b*sin(MyAngle); 
end

% Rotate and shift the samples
TrainingSamples=[cos(phi) -sin(phi);sin(phi) cos(phi)]*TrainingSamples+repmat(Center,[1 NumTrainSamples]);
TestSamples=[cos(phi) -sin(phi);sin(phi) cos(phi)]*TestSamples+repmat(Center,[1 NumTestSamples]);

% Collect the geometric and algebraic parameters
ParG=[Center(1) Center(2) a b phi]';
ParA=GtoA(ParG,1);

% Convert to natural form
ParN=GtoN(ParG);