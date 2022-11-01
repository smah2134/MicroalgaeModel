function [I] = LightModel_BeerLambert_Perpendicular_Tube(X,I_0,R,x_i,y_i)




k = 80;                                                  %% Absorption/Extinction coefficient

I = I_0.*exp(-k*X*((sqrt((R^2)-(x_i.^2))-y_i)/1000));    %% Light intensity at every point in the tube. X is the biomass concentration, sqrt((R^2 - x^2) y) is the mono-dimensional light path length in y-direction


end