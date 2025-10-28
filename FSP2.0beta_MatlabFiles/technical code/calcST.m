%Darren Pais
%Surveillance and Optimization, Reservoir
%2016+

% Outputs: S: Storativity
%          T: Transmissivity

% Inputs: h: reservoir thickness
%         phi: porosity fraction
%         kap: intrinsic permeability
%         valsarr: array of other inputs as below (to conform with rest of
%         code, a subset of these are used in the calc as below) 


function [S , T] = calcST(valsarr , h , phi , kap)

%all SI below with conversions as needed
rho = valsarr(5); %density
g = valsarr(6); %accl gravity
mu = valsarr(7); %dynamic viscosity
beta = valsarr(8); %fluid compressibility
alphav = valsarr(9); %vertical compressibility of aquifer

kap = kap*10^-3*9.9e-13; %convert permeability mD to m^2

S = rho*g*h*(alphav+phi*beta); %storativity 

% specific storage: (For MODFLOW Comparison
%  SpecificStorage = S./h  % rho*g*(alphav+phi*beta); % 

K = kap*rho*g/mu ;  %saturated hydraulic conductivity
T = K*h ;%transmissivity

end