function calldatHydroMc(~,~,hDV)

if hDV.plotdata.printFunctionName
    disp(['running calldatHydroMc'])
end

title = 'Uniform Distribution bounds' ;
% i=6 ; % where does the data go, see dataentry.m
noshow = [0] ; %what not to show

if ~isfield(hDV.data,'probHydrology')
    hDV.data.probHydrology={};
end

numVariables=8;
if isfield(hDV.data.probHydrology,'sigvals')
    
    vals = hDV.data.probHydrology.sigvals;

else
    vals=[zeros(numVariables-1,1);200];
    
%     vals = hDV.data.sigvals;  % backward compatible
%     hDV.data.distributions.vals = vals ;
%     
%     % 0=constant and noshow, 1=uniform distribution, 2=constant, ...
%     hDV.data.distributions.distriutionType=ones(12,1);
%     hDV.data.distributions.distriutionType(noshow,1)=0;
end


% % parameters that can vary
% % get storativity and transmissivity
% aqthick = hDV.data.reservoir.vals(1)*0.3048 ; %ft to meters
% porosityFraction=hDV.data.reservoir.vals(2)/100;
% perm_mD= hDV.data.reservoir.vals(3);
% rho = hDV.data.adv.vals(5)  ;% [kg/m^3]
% dynamicVisc=hDV.data.adv.vals(7) ;%[Pa.s]
% FluidCompressibility= hDV.data.adv.vals(8) ;%[Pa^-1]
% RockCompressibility= hDV.data.adv.vals(9) ;% [Pa^-1]

% this will go in dataentry
deterministicVals=[hDV.data.reservoir.vals(1);hDV.data.reservoir.vals(2);hDV.data.reservoir.vals(3);hDV.data.adv.vals(5);hDV.data.adv.vals(7);hDV.data.adv.vals(8);hDV.data.adv.vals(9);vals(8)];
txt2 = {'Aquifer Thickness [',' ft]' ;
    'Porosity [',' %]' ;
    'Perm [',' mD]' ;
    'fluid density [' ,' kg/(m^3)]';
    'dynamic viscosity [',' Pa.s]' ;
    'Fluid Compressibility [',' Pa^-1]' ;
    'Rock Compressibility [',' Pa^-1]';
    '#Hydrologic Iterations=',', change?'};

distirbutionTxt=cell(numVariables,1);
for kj=1:numVariables % concatenate string with value from input
    distirbutionTxt{kj,1}=[txt2{kj,1},num2str(deterministicVals(kj)),txt2{kj,2}];
end

% bound1OrPlusMinus2=[1;1;1;1;1;2;2;1;1;2;1;1]; % parameters that vary with each fault(2), or that vary with stresses(1)
% hDV.data.distributions.bound1OrPlusMinus2=bound1OrPlusMinus2;
dataentryProbabilisticHydrology(hDV,title,distirbutionTxt,vals,noshow);
hDV.data.probHydrology.distirbutionTxt=distirbutionTxt; % save text


end


