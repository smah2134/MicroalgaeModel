function [qSun] = solarRadiation(inputRadiation) 

%% Solar heat flux
%
%This function calculates the heat flux from solar radiation
%
%Assumptions: the measured global horizontal radiation hits the projected
%area of the tubes. cosine losses are neglected


cAbsorption = 0.95;             %Absorption coefficient of reactor water
cTransmission = 0.95;           % Transmission coeffcicient
surfaceP = 2.81793;             % Projected surface (m²)


qSun = inputRadiation * surfaceP * cTransmission * cAbsorption*10^-3;

end

