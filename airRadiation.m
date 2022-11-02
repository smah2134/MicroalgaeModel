function [qAir] = airRadiation(pVWater,AmbientTemperature,immersion)

%% Air radiation function
%
%This function calculates the thermal radiation from the surrounding air to
%the reactor surface



cBbRad = 5.67*10^-8;                                            % Stefan-Boltzmann constant
surfaceRad = 7.76;                                              % Outer total surface area (m²)
cAirRadiation = 1.048-10.^-((0.4* ((pVWater)).^0.2));



%% Air radiation based on the atmospheric or sky temperature
%
%Assumption: view factor for all immersion levels: 0.5
%Change view factors according to your caclulations/estimations

TSky = 0.0552.*(AmbientTemperature+273.15).^1.5;

if immersion == 0.5

    qAir = (cBbRad * ((TSky).^4) * (cAirRadiation) * surfaceRad*0.5)*10^-3;

elseif immersion == 0.25

    qAir = (cBbRad * ((TSky).^4) * (cAirRadiation) *surfaceRad*0.5)*10^-3;

else
    qAir = (cBbRad * ((TSky).^4) * (cAirRadiation) * surfaceRad*0.5)*10^-3;

end


%% Air radiation based on the water vapour pressure
%
%Alternative calculation method. View factors must be adjusted/estimated

% if immersion == 0.5
% 
%     qAir = (cBbRad * ((AmbientTemperature+273.15).^4) * (cAirRadiation) * surfaceRad*0.5*0.85)*10^-3;
% 
% elseif immersion == 0.25
% 
%     qAir = (cBbRad * ((AmbientTemperature+273.15).^4) * (cAirRadiation) *surfaceRad*0.667*0.6)*10^-3;
% 
% else
%     qAir = (cBbRad * ((AmbientTemperature+273.15).^4) * (cAirRadiation) * surfaceRad*0.5)*10^-3;
% 
% end



end

