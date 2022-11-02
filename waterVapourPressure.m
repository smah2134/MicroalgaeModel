function [wvpressure] = waterVapourPressure(AmbientTemperature)

%This function calculates the water vapour pressure for the calculation of
%the air radiation
%
%Only necessary when not using Tsky

%% Antoine Equation

A = 8.07131;
B = 1730.63;
C = 233.426;

wvpressure = (10.^(A-(B/(C + AmbientTemperature + 273.15)))).*133.322387415;


%% Simplified calculation

% wvpressure = (exp((20.386-(5132/(AmbientTemperature+273.15)))))*1.333224;



end