% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Pore pressure to failure Monte Carlo Calculation. Runs the mohrs_3D.m code
%with stochastic inputs

%inputs: indatacell: Cell array of baseline inputs into mohrs_3D.m
%        sigcell: Cell array of uniform standard deviations corresponding
%        to each input in indatacell
%        dpmax: Maximum dp range over which to run the calculation
%        nruns: How many runs for each set of inputs to compute results
function [outs,ins]  =ppfail_MC(indatacell,sigcell,hDV,nruns)


This function is not used now (and this line will lead to error if used) - SPL

if nargin < 4
    nruns = 1000 ; %default to 1000 runs if input not provided
end

[outs,ins] = monte_carlo(@mohrs_3D_v2,indatacell,sigcell,hDV,nruns) ; %performs MC calculation

end