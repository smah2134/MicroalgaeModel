%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       Main model for the pridiction of microalgal growth
%
%           Author: Alexander Hofmann
%
%           This model couples the light and temperature sub-models to the
%           kinetic growth sub-model to calculate the microalgae concentration
%           for eacht time step.
%       
%
%
%
%
%
%           Copyright 2022 Alexander Hofmann
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






clear
close all


%% Load Environmental Data

load 'filename'                      %% Change to file name (e.g. matlab table data format)
LightPAR= ;                          %% Specify PAR data (micromol/m²s)
LightGlobal= ;                       %% Specify global radiation data
AmbientTemperature= ;                %% Specify ambient temperature data
WaterTemperature = ;                 %% Specify water temperature data
ReactorTemperature= ;                %% Specify reactor temperature data
datetime =      ;                %% Date time     




%% Set initial Conditions



X_start = 0.5;                                                                       %Initial microalgae concentration (g/L)
t_total = 11;                                                                        %Total simulation time (d) 
timestep = 0.5/(24*60);                                                              %Time step size (d)
timestep_s =timestep*24*3600;                                                    %Time step size (s)               
timesteps_total = t_total/timestep;                                                  %Total time steps                                                                                
X(1) = X_start_AH1;                                                                  %Microalgae concentration at t=1
    

immersion = 0.5;                                                                    %Immersion level, possible values: 0.5, 0.25, 0
TReactorOld = ReactorTemperature(1);                                                %Start temperature for reactor temperature simulation


%% Load geometry

[R,V,x_i,y_i,V_total] = ReactorGeometry_Tube;                                                    %Load tube geometry function



for t=1:1:timesteps_total


    %% Read Light Data

    %Read incident light intensity data for specific time step
    I_0(t) = LightPAR(t,1);                                  %Incident light (PAR) intensity (micromol/m²s)



    %% Read Temperature Data
    %Read reactor temperature data for specific stime step
    T_R(t) = ReactorTemperature(t,1);   %Measured reactor temperature (°C)


    %Reactor temperature simulation

    [T_Rsim(t+1), TReactorOld(t+1)] = ReactorTempSimMain(immersion,LightGlobal(t),AmbientTemperature(t),WaterTemperature(t),TReactorOld(t));

    %% Light Distribution Model (comment/uncomment for use)
        
       % Beer-Lambert Model
        [I]=LightModel_BeerLambert_Perpendicular_Tube(X(t),I_0(t),R,x_i,y_i);  %%Light intensity at every grid point in micromol/m²s



    %% Growth Kinetics Model (comment/uncomment for use)


        [mu_max,mu,X_dv] = GrowthModel_simpleMonod_CTMI_Resp(T_R(t),I,I_0(t),timestep,X(t));  %%Monod model combined with cardinal temperature model with inflexion. change T_R to T_Rsim for the simulated temperature



        Biomass_abs = X_dv.*V;                                  % New biomass of each volume element
        X(t+1) =  sum(Biomass_abs,'all')./V_total;              % New biomass concentration of each volume element
        time(t) = t.*timestep;                                  % Real time
        timeO(t+1) = t.*timestep;                               % Real time with off-set


end


%% Plot results

sz = 10;

figure(1)
plot(timeO,X, 'black')
legend({'Simulation'},'Location','northwest','Orientation','vertical')
xlabel('Cultivation time (d)')
ylabel('Biomass concentration (g/L)')
