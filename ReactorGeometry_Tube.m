function [R,V,x_i,y_i,V_total] = ReactorGeometry_Tube

%% Geometry function
%
%This function creates a 2d mesh for the calculation of the light
%distribution within the pipes


%% Geometry dimensions

R = 30.4;                           %Inner tube radius (mm)
dx = 0.1;                           %Spatial discretization - radius
dtheta = 1;                         %Spatial discretization - angle (°)
dthetaR = dtheta*(pi/180);          %Conversion degree - radiant
L = 54.9284*1000;                   %Total tube lenght (mm) 54.9284
x = ((0+(dx/2)):dx:(30.4-(dx/2)));  %Interval for mesh generation. Grid points are assumed to be the mid-points of differential volume elements
rho = x;
dx1 = -dx;
Rv = (0:dx1:30.4);
Rstep = R/dx;

%%Mesh generation

theta = ((0+(dtheta/2)):dtheta:(360-(dtheta/2)))*pi/180;
thetasteps = 360/dtheta;
[th, r] = meshgrid(theta, rho);
y_i = r.*sin(theta);                %Conversion of cylindrical coordinates to Cartesian
x_i = r.*cos(theta);

x = (dx:dx:(30.4));  %Interval for mesh generation
rho = x;
theta = (dtheta:dtheta:(360))*pi/180;
[th, r] = meshgrid(theta, rho);

 V = r.*dx.*dthetaR*L.*10^-9;
 V_total = sum(V,"all");


end