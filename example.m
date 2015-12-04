clear all
close all

dir = 'H:\XRayAA_WaveLets\myCustomMatlabFunction\xraylibForGitHub\readyForGitHub';
cd(dir)

mex -IC:\boost_1_57_0 mexCinCoeffThinNoAbs.cpp Element.cpp XRLFunctions.cpp ElementEnergyLine.cpp Compound.cpp IncidentEnergyObj.cpp  CinFunctioncs.cpp
mex -IC:\boost_1_57_0 mexCinCoeff.cpp Element.cpp XRLFunctions.cpp ElementEnergyLine.cpp Compound.cpp IncidentEnergyObj.cpp  CinFunctioncs.cpp

zA = [29, 79];
wA = [0.4915166904661467, 0.5084833095338532];
lambdaLines = [0, 2];

%{
mZ1 = 72.58999633789062; % Ge
mZ2 = 63.540000915527344; % Cu
x = 0.02;
massSum = (1-x)*mZ1 + x*mZ2;
WZ1 = (1-x)*mZ1/massSum;
WZ2 = x*mZ2/massSum;

zA = [32, 29];
wA = [WZ1 WZ2];
%}

N = length(zA);

rho = 11.5; % g/cm^3
mT = 200*rho/10000;

theta2 = 0*pi/180;
theta1 = 45*pi/180;
energy = 17.44;

CThin = zeros(N,N);
C = zeros(N,N);

for aa = 1:N
    for bb = 1:N
        CThin(aa, bb) = mexCinCoeffThinNoAbs([zA(aa) wA(aa)],[zA(bb) wA(bb)], energy, lambdaLines(aa), zA, wA, theta1, theta2, mT);
        C(aa, bb) = mexCinCoeff([zA(aa) wA(aa)],[zA(bb) wA(bb)], energy, lambdaLines(aa), zA, wA, theta1, theta2);
    end
end

% thin C matrix
CThin

% thick C matrix
C

