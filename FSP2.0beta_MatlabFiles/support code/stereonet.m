% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface


%plot stereonet projection for faults given dip and strike
%output: Theta: Polar coordiate angle of stereonet curve
%        Rho: Polar coordiante distance (rho) of stereonet curve
%        poleth: Polar corrdinate angle for pole
%        polerho: Polar coordinate radius for pole

function [Theta, Rho, poleth, polerho]=stereonet(dip,strike)

if nargin==0
    dip=    [45 ; 90 ; 10] ;
    strike =[45 ; 60 ; 90] ;   
end

strike = strike*pi/180 ; %to radians
dip    = dip*pi/180 ; %to radians

%define pole
poleth = strike+pi/2;
polerho = dip/(pi/2) ; 


R=1 ;
rake=0:pi/180:pi;
% rake=pi:pi/180:2*pi;

%placeholder
Theta=zeros(length(dip),length(rake)) ; Rho=Theta ;

%calculate projections
for i=1:1:length(strike)
    plunge = asin(sin(dip(i)).*sin(rake));
    trend  = strike(i) + atan2(cos(dip(i)).*sin(rake),cos(rake)) ;
    rho=R.*tan(pi/4 - plunge/2) ;
%     Theta(i,:) = trend;  %as per MATLAB goes CCW from 0300 .. to go CW from 1200 use pi/2-trend
    Theta(i,:) = trend+pi;  %make lower hemisphere not upper hemisphere? -FRW
    Rho(i,:) = rho ;
end

%plotting
if nargin==0
    for i=1:1:length(strike)
          hp=polar(Theta(i,:),Rho(i,:)) ;

        set(hp,'linewidth',2);
        hold all;
    end
end

end