function [qTransfer,qTransferWater,qTransferAir,alpha_pipe_outer_water,NU_pipe_water_free,Ra_pipe_water,L_Ra_pipe_water] = heatTransferAirAndWater(immersion,tReactor, tOutside, tWater,b)  


%% Heat transfer between the microalgae culture and the surrounding air and water
%
%This function calculates the heat flux from/to the microalgae culture
%from/to the surrounding air and water
%
%
%A constant air and water flow speed is assumed

A_o = 7.76;                         % Outer reactor surface (m²)
A_i = 7.49;                         % Inner reactor surface (m²)
cHeatTransferPVC = 0.15;            % Thermal conductivity PVC (W/mK)
%A_m_pipe = 7.986;
Re_pipe = 42560;                    % @ flow velocity=0.7 m/s
Pr_WaterPipe = 6.95;
d_i = 60.8*10^-3;
d_o = 63*10^-3;
L_pipe = 41.06;

hCond_water = 0.6;
hCond_air = 0.0262;
s = 0.0011;
Visc_water_dyn = 10^(-6);
Visc_air_dyn = 18*10^(-6);
thermalDiff_water = hCond_water/(997*4184);
thermalDiff_air = hCond_air/(1.225*1005);
L_Ra_pipe_water = 0.5*d_o*(pi);
g = 9.81;
betaWater = 0.1455*10^-3;
betaAir = 1/(tOutside+273.15);
Pr_pipe_free_water = Visc_water_dyn/thermalDiff_water;
fPr = (1+((0.559/Pr_pipe_free_water)^(9/16)))^(-(16/9));




Pr_pipe_free_air = Visc_air_dyn/thermalDiff_air;
fPr_air = (1+((0.559/Pr_pipe_free_air)^(9/16)))^(-(16/9));
L_Pipe_water_forced_outer = 0.5*pi*d_o;
L_Ra_pipe_water_25 = (1/3)*pi*d_o;
L_Ra_pipe_air_25 = (2/3)*pi*d_o;
L_Ra_pipe_air_0 = 0.5*pi*d_o;


A_o_05_water = A_o*0.5;
A_o_05_air = A_o_05_water;
A_o_025_water = 0.33*A_o;
A_o_025_air = A_o - A_o_025_water;

A_i_05_water = A_i*0.5;
A_i_05_air = A_i_05_water;
A_i_025_water = A_i*0.33;
A_i_025_air = A_i - A_i_025_water;

    

Xi_pipe = (1.8*(log10(Re_pipe)-1.5))^-2;
NU_pipe = (((Xi_pipe/8)*Re_pipe*Pr_WaterPipe)/(1+(12.7*sqrt(Xi_pipe)*((Pr_WaterPipe^(2/3))-1))))*(1+((d_i/L_pipe)^(2/3)));
alpha_pipe_inner = (NU_pipe*hCond_water)/d_i;


v_water_outer = 0.2;
v_air_outer = 2;




if immersion == 0.5

    %% Dimensionless numbers and geometry

    %% Pipe flow 
  

    %% Heat transfer water

     L_char_w = 0.5*pi*d_o;

     % Forced outer convection
     a_pipes = 270/63;
     %Psi_pipes = 1-(pi/(4*a_pipes));
     Re_pipe_forced_outer_water = (v_water_outer*L_char_w)/(Visc_water_dyn);
     NU_pipe_outer_lam_water = 0.664*(sqrt(Re_pipe_forced_outer_water))*(Pr_pipe_free_water^(1/3));
     NU_pipe_outer_tub_water = (0.037*(Re_pipe_forced_outer_water^(0.8))*Pr_pipe_free_water)/(1 + (2.443*(Re_pipe_forced_outer_water^(-0.1)))*((Pr_pipe_free_water^(2/3))-1));
     NU_pipe_water_forced_outer = 0.3 + sqrt((NU_pipe_outer_lam_water^(2))+(NU_pipe_outer_tub_water^(2)));%850;
     
     % Free outer convection
     Ra_pipe_water = (betaWater.*g.*(abs((tReactor-tWater))).*(L_Ra_pipe_water.^3))./(Visc_water_dyn.*thermalDiff_water);
     NU_pipe_water_free = (0.752 + (0.387.*((Ra_pipe_water*fPr).^(1/6)))).^2;
     NU_pipe_water_mixed = ((NU_pipe_water_forced_outer^3) + (NU_pipe_water_free^3))^(1/3);
     alpha_pipe_outer_water = (NU_pipe_water_mixed*hCond_water)/L_char_w;
     

    %% Heat transfer air

    Ra_pipe_air = (betaAir.*g.*(abs((tReactor-tOutside))).*(L_Ra_pipe_water.^3))./(Visc_air_dyn.*thermalDiff_air);
    NU_pipe_air_free = (0.752 + (0.387*((Ra_pipe_air*fPr_air)^(1/6))))^2;
    Re_pipe_forced_outer_air = (v_air_outer*L_char_w)/(Visc_air_dyn);

    NU_pipe_outer_lam_air = 0.664*(sqrt(Re_pipe_forced_outer_air))*(Pr_pipe_free_air^(1/3));
    NU_pipe_outer_tub_air = (0.037*(Re_pipe_forced_outer_air^(0.8))*Pr_pipe_free_air)/(1 + (2.443*(Re_pipe_forced_outer_air^(-0.1)))*((Pr_pipe_free_air^(2/3))-1));
    NU_pipe_air_forced_outer = 0.3 + sqrt((NU_pipe_outer_lam_air^(2))+(NU_pipe_outer_tub_air^(2)));%850;
    NU_pipe_air_mixed = ((NU_pipe_air_forced_outer^3) + (NU_pipe_air_free^3))^(1/3);
  
    alpha_pipe_outer_air = (NU_pipe_air_mixed*hCond_air)/L_Ra_pipe_water; %175;


    %% Total heat transfer

    A_m_pipe = (A_o_05_air - A_i_05_water)/log(A_o_05_water/A_i_05_water);

    k_w = (1/A_o_05_water)*(1/((1/(alpha_pipe_outer_water*A_o_05_water))+(s/(cHeatTransferPVC*A_m_pipe))+(1/(alpha_pipe_inner*A_i_05_water))));

    k_a = (1/A_o_05_air)*(1/((1/(alpha_pipe_outer_air*A_o_05_air))+(s/(cHeatTransferPVC*A_m_pipe))+(1/(alpha_pipe_inner*A_i_05_water))));

    qTransferWater = k_w*A_o_05_water*(tWater-tReactor);

    qTransferAir = k_a*A_o_05_air*(tOutside - tReactor);




elseif immersion == 0.25

    %% Dimensionless numbers and geometry

%     %% Pipe flow 
%     L_char_w = 
% 
%     %% Heat transfer water
% 
     L_char_w = (1/3)*0.5*pi*d_o;
     a_pipes = 270/63;
     Psi_pipes = 1-(pi/(4*a_pipes));
     Re_pipe_forced_outer_water = (v_water_outer*L_char_w)/(Psi_pipes*Visc_water_dyn);
     NU_pipe_outer_lam_water = 0.664*(sqrt(Re_pipe_forced_outer_water))*(Pr_pipe_free_water^(1/3));
     NU_pipe_outer_tub_water = (0.037*(Re_pipe_forced_outer_water^(0.8))*Pr_pipe_free_water)/(1 + (2.443*(Re_pipe_forced_outer_water^(-0.1)))*((Pr_pipe_free_water^(2/3))-1));
     NU_pipe_water_forced_outer = 0.3 + sqrt((NU_pipe_outer_lam_water^(2))+(NU_pipe_outer_tub_water^(2)));%850;

     Ra_pipe_water = (betaWater.*g.*(abs((tReactor-tWater))).*(L_Ra_pipe_water_25.^3))./(Visc_water_dyn.*thermalDiff_water);
     NU_pipe_water_free = (0.752 + (0.387*((Ra_pipe_water*fPr)^(1/6))))^2;
     NU_pipe_water_mixed = ((NU_pipe_water_forced_outer^3) + (NU_pipe_water_free^3))^(1/3);
     alpha_pipe_outer_water = (NU_pipe_water_mixed*hCond_water)/L_char_w;

   

    %% Heat transfer air

    Ra_pipe_air = (betaAir.*g.*(abs((tReactor-tOutside))).*(L_Ra_pipe_air_25.^3))./(Visc_air_dyn.*thermalDiff_air);
    NU_pipe_air_free = (0.752 + (0.387*((Ra_pipe_air*fPr_air)^(1/6))))^2;
    Re_pipe_forced_outer_air = (v_air_outer*L_char_w)/(Visc_air_dyn);

    NU_pipe_outer_lam_air = 0.664*(sqrt(Re_pipe_forced_outer_air))*(Pr_pipe_free_air^(1/3));
    NU_pipe_outer_tub_air = (0.037*(Re_pipe_forced_outer_air^(0.8))*Pr_pipe_free_air)/(1 + (2.443*(Re_pipe_forced_outer_air^(-0.1)))*((Pr_pipe_free_air^(2/3))-1));
    NU_pipe_air_forced_outer = 0.3 + sqrt((NU_pipe_outer_lam_air^(2))+(NU_pipe_outer_tub_air^(2)));%850;
    NU_pipe_air_mixed = ((NU_pipe_air_forced_outer^3) + (NU_pipe_air_free^3))^(1/3);
  
    alpha_pipe_outer_air = (NU_pipe_air_mixed*hCond_air)/L_Ra_pipe_air_25; %175;
     

    A_m_pipe = (A_o_025_water - A_i_025_water)/log(A_o_025_water/A_i_025_water);

    A_m_pipe_air = (A_o_025_air - A_i_025_air)/log(A_o_025_air/A_i_025_air);

    k_w = (1/A_o_025_water)*(1/((1/(alpha_pipe_outer_water*A_o_025_water))+(s/(cHeatTransferPVC*A_m_pipe))+(1/(alpha_pipe_inner*A_i_025_water))));

    k_a = (1/A_o_025_air)*(1/((1/(alpha_pipe_outer_air*A_o_025_air))+(s/(cHeatTransferPVC*A_m_pipe_air))+(1/(alpha_pipe_inner*A_i_025_water))));

    qTransferWater = k_w*A_o_025_water*(tWater-tReactor);

    qTransferAir = k_a*A_o_025_air*(tOutside - tReactor);

else


    Ra_pipe_air = (betaAir.*g.*(abs((tReactor-tOutside))).*(L_Ra_pipe_air_0.^3))./(Visc_air_dyn.*thermalDiff_air);
    NU_pipe_air_free = (0.752 + (0.387*((Ra_pipe_air*fPr_air)^(1/6))))^2;
    Re_pipe_forced_outer_air = (v_air_outer*L_Ra_pipe_air_0)/(Visc_air_dyn);

    NU_pipe_outer_lam_air = 0.664*(sqrt(Re_pipe_forced_outer_air))*(Pr_pipe_free_air^(1/3));
    NU_pipe_outer_tub_air = (0.037*(Re_pipe_forced_outer_air^(0.8))*Pr_pipe_free_air)/(1 + (2.443*(Re_pipe_forced_outer_air^(-0.1)))*((Pr_pipe_free_air^(2/3))-1));
    NU_pipe_air_forced_outer = 0.3 + sqrt((NU_pipe_outer_lam_air^(2))+(NU_pipe_outer_tub_air^(2)));%850;
    NU_pipe_air_mixed = ((NU_pipe_air_forced_outer^3) + (NU_pipe_air_free^3))^(1/3);
  
    alpha_pipe_outer_air = (NU_pipe_air_mixed*hCond_air)/L_Ra_pipe_air_0; 

    A_m_pipe = (A_o - A_i)/log(A_o/A_i);

    k_a = (1/A_o)*(1/((1/(alpha_pipe_outer_air*A_o))+(s/(cHeatTransferPVC*A_m_pipe))+(1/(alpha_pipe_inner*A_i))));

    qTransferWater = 0;

    qTransferAir = k_a*A_o*(tOutside - tReactor);

end





qTransfer =  (qTransferAir+qTransferWater)*10^-3;
end

