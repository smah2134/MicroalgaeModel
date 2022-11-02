function [qWater] = waterRadiation(TWater,immersion)

%% Radiation from the water body
%
%This function estimates the heat flux radiataed from the water body
%surrounding the tubes
%
%Assumption: thermal radiation from the water body acts on the immersed
%part of the tube and on the emerged tube part from its surface
%
%View factors are estimated

cBbRad = 5.67*(10^(-8));         % [W/(m2K4] Stefan-Boltzmann Konstante für Radiation eines Black Body (SIGMA)
cTransmission =0.95*10^-3;       % Transmission coeffcicient
emissivityW = 0.95;              % Emissivity water
surfaceR = 7.76;                 % Outer reactor surface (m²)


if immersion == 0.5

qWaterImmersed = cBbRad * (emissivityW.^2) * cTransmission * ((TWater+273.15).^4) * surfaceR*0.5;  

qWaterAir = cBbRad * (emissivityW.^2) * cTransmission * ((TWater+273.15).^4) * surfaceR*0.5*0.1;    

elseif immersion == 0.25

qWaterImmersed = cBbRad * (emissivityW.^2) * cTransmission * ((TWater+273.15).^4) * (1/3)*surfaceR;  

qWaterAir = cBbRad * (emissivityW.^2) * cTransmission * ((TWater+273.15).^4) * (2/3)*surfaceR*0.2;    

else 

qWaterImmersed = 0;       
qWaterAir = cBbRad * (emissivityW.^2) * cTransmission * ((TWater+273.15).^4) * surfaceR*0.6;   %estimated view factor to account for reflections

end



%qWater =  immersion *qWaterEmerged + (1-immersion) * qWaterAir;
qWater =  qWaterImmersed + qWaterAir;

end

