function [TReactorSim, TReactorOld] = ReactorTempSimMain(immersion,LightGlobal,AmbientTemperature,WaterTemperature,TReactorOld)

%% This function calculates the reactor temperature depending on the immersion level, the global radiation intensity in W/m², the ambient temperature, and the air temperature
%%
%% This model was created by Pia Leminski in 2022 and later modified by Alexander Hofmann

%Assumptions: heat flows from friction and microalgae metabolism are
%neglected

% Assign start temperature
 %TReactorOld(1) = Tdata(1);
TReactorSim = (TReactorOld);
% 
Vculture = 0.16; % (m3)




heatcapacitywater = 4.184;
waterdensity = 997;

    qTot =(solarRadiation(LightGlobal)) +(reactorRadiation(TReactorOld)) +(airRadiation(waterVapourPressure(AmbientTemperature),AmbientTemperature,immersion)) + (waterRadiation(WaterTemperature, immersion)) +(heatTransferAirAndWater(immersion, TReactorOld, AmbientTemperature, WaterTemperature));
    TReactorNew = TReactorOld+ (qTot*30 / (Vculture*heatcapacitywater*waterdensity));
     
    TReactorOld = TReactorNew;
    TReactorSim= TReactorOld;

 end