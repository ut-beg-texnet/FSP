% Rall Walsh
% run probabilistic hydrology
% runs the hydro1D code
%with stochastic inputs

%inputs: indatacell: Cell array of baseline inputs into hydro1D.m
%        sigcell: Cell array of high-low bounds
%        nruns: How many runs for each set of inputs to compute results?

% output
%        1 column per fault
% value is Pore pressure on that fault in PSI

function [outs,ins] =ppOnFault_MC(indatacell,sigcell,nruns)

if nargin == 0
    Xwell=[0.4;10.3];
    Ywell=[0;0];
    % constant well rate data
    allWellsDatenumBarrelsPerDay{1,1}=[   737061.000694444  ,     0;
        749479.000694444  ,       23000];
    allWellsDatenumBarrelsPerDay{2,1}=[   737061.000694444  ,     0;
        749479.000694444  ,       12500];
    xFault=[4,18,7,24]';yFault=zeros(size(xFault));
    tSlider = 2025; % get(hDV.hdsldr(1),'value') ; % time
    wells = 1:2; %indices of wells to use
    g =9.81  ; % constant [m/s^2]
    
    % parameters that can vary
    % get storativity and transmissivity
    aqthick = 100*0.3048 ; %100 ft to 30.48 meters
    porosityFraction=0.1;
    perm_mD= 200;
    rho = 1000 ; %fluid density [kg/m^3]
    dynamicVisc= 0.0008 ;%[Pa.s]
    FluidCompressibility=  3.6e-10 ;%[Pa^-1]
    RockCompressibility=  1.08e-09 ;% [Pa^-1]
    
    indatacell={Xwell,Ywell,wells,xFault,yFault,g,tSlider,aqthick,porosityFraction,perm_mD,rho,...
        dynamicVisc,FluidCompressibility,RockCompressibility,{allWellsDatenumBarrelsPerDay}};
    % vary perm and porosity
    sigcell ={0,0,0,0,0,0,0,0,.05,150,0,0,0,0,0} ;
    nruns =200 ;
elseif nargin == 2
    nruns = 200 ; %default to 1000 runs if input not provided
end

% each column is a fault in Outs
[outs,ins] = monte_carloHydro(@hydro1D,indatacell,sigcell,nruns,length(indatacell{4}));  %performs MC calculation

if nargin==0
    figure ;
    hold on
    %     outs(end-500:end)=max(outs) % test distribution reflected in CDF
    for kjj4=1:size(outs,2)
        [f,x] = ecdf(outs(:,kjj4)) ;
        plot(x,1-f,'linewidth',2) ; % probability of not exceeding
        xlabel('PP on fault') ;
        ylabel('CDF Probability of not exceeding [PSI]') ;
    end
    grid on
    set(gcf,'color','w')

end

end