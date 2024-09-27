% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface
% Setup data input
% modified by Rall
% Modified by Suvrat Lele

function setupdata(hDV)

hDV.data.nwells = 1;

hDV.data.stress = [];
hDV.data.reservoir = [];
hDV.data.inject = [];
hDV.data.fault = [];
hDV.data.realWellData = [];
hDV.data.probHydrology =[];
hDV.plotdata.hydprob.yearDataShown=NaN;

hDV.data.stress.txt = {
    'Vertical Stress Gradient [psi/ft]' ;...
    'Min Horiz. Stress Gradient [psi/ft]'  ;...
    'Max Horiz. Stress Gradient [psi/ft]'  ;...
    'Max Hor Stress Direction [deg N CW]' ; ...
    'Reference Depth for Calculations [ft]'; ...
    'Initial Res. Pressure Gradient [psi/ft]'};
hDV.data.stress.vals = zeros(length(hDV.data.stress.txt),1) ;
hDV.data.stress.aphi.use = 0;
hDV.data.stress.aphi.vals = [0 0] ;
hDV.data.stress.aphi.sigvals = [0,0];

hDV.data.reservoir.txt = { 'Aquifer Thickness [ft]'  ;...
    'Porosity [%]'  ; ...
    'Permeability [mD]';...
    'Stochastic Mag [psi]';...
    'Poisson''s Ratio'
    } ;
hDV.data.reservoir.vals = zeros(length(hDV.data.reservoir.txt),1) ;
hDV.data.reservoir.vals(5) = 0.5 ; % hard code the poisson's ratio for now
hDV.data.reservoir.loadedHeaderLines=1;
hDV.data.reservoir.importHydrology=0; % 0 does hydrology calculation internally, 1 imports XYZPT model like from MODFLOW

hDV.data.inject.txt = {'x [km]' ;...
    'y [km]'  ;...
    'Inj. Rate [bbl/day]'  ;...
    'Start Year [yr]' ;...
    'End Year [yr]'} ;
hDV.data.inject.vals = zeros(hDV.data.nwells_max,length(hDV.data.inject.txt)) ;

hDV.data.realWellData.inputStringColumnNames={'UniqueID/Name';'Easting (km)';'Northing (km)';'Year';'Month (1-12)';'InjectionVolume (bbl/month)'};
hDV.data.realWellData.stringsWellDataAdvanced=cell(2,size(hDV.data.realWellData.inputStringColumnNames,1));
hDV.data.realWellData.stringsWellDataAdvanced=[hDV.data.realWellData.stringsWellDataAdvanced{:,:}];

% hDV.data.realWellData.inputStringVariableNames={'API';'XEast';'YNorth';'Year';'monthNumber';'InjectionVolume';'Pressure';'WellName';'WellType'};
hDV.data.realWellData.columnIsNumber=logical([0,1,1,1,1,1]); % 1 for columns that take numeric data only, 0 for columns that can take letters and numbers
hDV.data.realWellData.loadedHeaderLines=1;



hDV.data.fault.txt = {['Number of faults (max ',num2str(hDV.data.NFAULTSMAX),')'] ;...
    'Friction Coefficient mu';...
    'Strike minimum [deg]';...
    'Strike maximum [deg]';...
    'Dip minimum [deg]';...
    'Dip maximum [deg]';}  ;
hDV.data.fault.vals = [20 0.6 0 0 0 0]';


hDV.data.adv.txt ={'Min x [km]' ; 'Max x [km]' ; 'Min y [km]' ; 'Max y [km]' ;...
    'Density [kg/m^3]' ; 'Accl Gravity [m/s^2]' ; 'Dynamic Viscosity [Pa.s]' ;...
    'Fluid Compressibility [Pa^-1]' ; 'Rock Compressibility [Pa^-1]' ;'Set Random Seed?' } ;
hDV.data.adv.vals = [0 20 0 20 1000 9.80665 .8*10^-3 3.6e-10   3*3.6e-10 0]'  ;  %all SI units

hDV.plotdata.minint=0;
hDV.plotdata.maxint=1;

if isdeployed 
    hDV.plotdata.printFunctionName=0;
else % print some function names to screen
     hDV.plotdata.printFunctionName=1;
end

end