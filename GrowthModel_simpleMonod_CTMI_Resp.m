function [mu_max,mu,X_dv,TempFactor] = GrowthModel_simpleMonod_CTMI_Resp(T_R,I,I_0,timestep,X)

%% Growth model function
%This function calculates the light and temperature-dependent specific growth rate and the microalgae concentration for every point in the tube 
%Values used here were determined for Tetradesmus obliquus
%
%
%Change values according to species characteristics

%% Cardinal Temperature with Inflexion

 %% Temperature dependent growth rate (based on Monod model)
 mu_max = 2.885;                                                    %(d^-1)
 Tmin = 0.65;                                                       %(°C)
 Tmax = 40.2;
 Topt = 37.64;
 TempFactor = ((T_R-Tmax).*(T_R-Tmin).^2)./((Topt-Tmin).*(((Topt-Tmin).*(T_R-Topt))-((Topt-Tmax).*(Topt+Tmin-(2.*T_R)))));

 %% Temperature dependent half-saturation constant

 K_s_max = 942.3;                                                   %micromol/m²s

 TmaxKs = 41.54;
 TminKs = 14.34;
 ToptKs = 35.06;
 TempFactorKs = ((T_R-TmaxKs).*(T_R-TminKs).^2)./((ToptKs-TminKs).*(((ToptKs-TminKs).*(T_R-ToptKs))-((ToptKs-TmaxKs).*(ToptKs+TminKs-(2.*T_R)))));



%% Respiration rate during the dark phase

%Assumption: no resperiation, because net growth rate data is used
%Change according to species characteristics

if I_0 > 50

          mu_dark = 0;
    
else
  
         mu_dark = -0.00;                                           %(d^-1)
end


%% Calculate half saturation constant

if TempFactorKs >0
    Ks = K_s_max.*TempFactorKs;

else
    Ks=0;

end

%% Calculate maximum growth rate

    mumax = mu_max.*TempFactor;


%% Calculate growth rate of every volume element
                                                  
    mu = (((mumax.*I)./(Ks+I)))+mu_dark;

%% Microalgae concentration for each volume element at timestep t
    
    X_dv = X.*exp(mu.*timestep);                                                        %New concentration of each volume element



end