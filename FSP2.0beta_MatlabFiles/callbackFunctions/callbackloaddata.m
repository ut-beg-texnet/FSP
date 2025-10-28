% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

% Load default data button callback

function callbackloaddata(~,~,hDV)

% make sure these match what's in surfaceviz
hDV.data.nwells_max = hDV.plotdata.wellsMaxFaultsMax(1) ; 
hDV.data.NFAULTSMAX = hDV.plotdata.wellsMaxFaultsMax(2);

            
hDV.data.stress.vals = [1.1 0.693 1.22 70 11000 0.43]' ;
hDV.data.reservoir.vals = [100 10 200 0*14.5 0.5]' ;
hDV.data.fault.vals = [20 .58 240 330 45 90]' ; hDV.data.fault.intype=1; hDV.data.fault.file=NaN;
hDV.data.nwells = 4 ;
hDV.data.inject.vals(1,:) = [17 12 27000  2015 2022];
hDV.data.inject.vals(2,:) = [3  4 23000  2016 2026] ;
hDV.data.inject.vals(3,:) = [14  7 15000  2017 2027] ;
hDV.data.inject.vals(4,:) = [9  9 12000  2018 2025] ;

% hDV.data.nwells = 2 ;
% hDV.data.inject.vals(1,:) = [16 12 270000  2015 2018];
% hDV.data.inject.vals(2,:) = [3  5 230000  2015 2017] ;


hDV.data.realWellData.use =  0 ;
hDV.data.adv.vals = [ -5 24 -3 20 1000 9.80665 .8*10^-3 3.6e-10   3*3.6e-10   20]' ; %region and fluid properties
% 2 and 4, 12 17 19 re good seeds?
% 20 good seed
hDV.data.probHydrology.sigvals=[25;3;50;2;0.000001;1e-12;1e-11;750];
hDV.data.probHydrology.probabilistic1VsDeterministic2=1;
hDV.data.probHydrology.distirbutionTxt={ 'Aquifer Thickness [100 ft]';...
    'Porosity [10 %]';...
    'Perm [200 mD]';...
    'fluid density [1000 kg/(m^3)]';...
    'dynamic viscosity [0.0008 Pa.s]';...
    'Fluid Compressibility [3.6e-10 Pa^-1]';...
    'Rock Compressibility [1.08e-09 Pa^-1]'};
% make datenumBarrelsPerDay for constant rate wells.
% YearEndOfModel=get(hDV.hdsldr(1),'max');% end of scaleruler.
% desired format for constant rate wells: thisWellDatenumBarrelsPerDay =
% [ start_datenumber_of_injection,   0; ...
%    end_datenumber_of_injection,   well_rate_bbls_per_day] ;
for kk=1:hDV.data.nwells
    hDV.data.inject.datenumBarrelsPerDay{kk,1}=[datenum(hDV.data.inject.vals(kk,4),1,1,0,1,0), 0; datenum(hDV.data.inject.vals(kk,5),1,1,0,1,0),hDV.data.inject.vals(kk,3)] ;
    
end
hDV.data.inject.datenumBarrelsPerDay=hDV.data.inject.datenumBarrelsPerDay(1:hDV.data.nwells,1);    % clear any old wells out



%sigma model for monte carlo default
hDV.data.sigvals = [.01,.01,.01,0,.01,5,5,5,0,.01,0,0]' ;

hDV.data.stress.aphi.use = 0;
hDV.data.stress.aphi.vals = [1.2277 0]';

% probabilistic trafficlight
hDV.plotdata.minint=0;
hDV.plotdata.maxint=0.5;
set(hDV.plotdata.pint.cmintxt,'string','0');
set(hDV.plotdata.pint.cmaxtxt,'string','.5');


% set(hDV.plotdata.pffot.cmintxt,'string',num2str(0));
% set(hDV.plotdata.pffot.cmaxtxt,'string',num2str(4000));

    
%load faults from random distribution
nfaults = hDV.data.fault.vals(1); %number of faults
Xmin=hDV.data.adv.vals(1); Xmax=hDV.data.adv.vals(2);
Ymin=hDV.data.adv.vals(3); Ymax=hDV.data.adv.vals(4);
Tmin=hDV.data.fault.vals(3);  Tmax=hDV.data.fault.vals(4);
Dmin=hDV.data.fault.vals(5);  Dmax=hDV.data.fault.vals(6);

hDV.data.fault.xf=Xmin+(Xmax-Xmin)*rand(nfaults,1);
hDV.data.fault.yf=Ymin+(Ymax-Ymin)*rand(nfaults,1);
hDV.data.fault.thf=Tmin+(Tmax-Tmin)*rand(nfaults,1); %fault orientation
hDV.data.fault.dipf =Dmin+(Dmax-Dmin)*rand(nfaults,1); %fault dips
hDV.data.fault.lenf=2+zeros(nfaults,1) ;
hDV.data.fault.muf=hDV.data.fault.vals(2) + zeros(nfaults,1) ;

set(hDV.bCalc,'enable','on') ;
callbackcalc([],[],hDV) % run calculation
end