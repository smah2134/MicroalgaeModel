function [qReactor] = reactorRadiation(TReactor) 

%% Reactor radiation
%
%This function calculates the thermal radiation losses of the reactor
% Calculates Reactor radiation 


cTransmission = 0.95*10^-3;                 % Transmission coeffcicient
cBbRad = 5.67*(10^(-8));                    % [W/(m2K4] Stefan-Boltzmann Konstante für Radiation eines Black Body (SIGMA)
surfaceRad = 7.76;%7.49;                    % whole inner surface
emissivityW = 0.95;                         % Emissivity water

qReactor = -((TReactor+273.15).^4) * cTransmission * cBbRad * surfaceRad * emissivityW;

end

