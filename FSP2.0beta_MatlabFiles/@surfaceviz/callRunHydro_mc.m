% Callback for running monte carlo hydrology
% modified by Rall Walsh from callrunmc.m by Darren Pais, ExxonMobil Upstream Research Company 2016
% Stanford


function callRunHydro_mc(~,~,hDV)

if nargin==0
    hDV=evalin('base', 'hSV'); % get hDV to test?
end

if hDV.plotdata.printFunctionName
    disp(['running callRunHydro_mc '])
end


[flag] = checkProbabilisticEntriesHyd(hDV); % chech hydrologic uncertainties and data again
if ~flag % don't calculate if data issue
    return
end
%
% % if calculate button red, run that first
% if all(get(hDV.bCalc,'backgroundcolor')==hDV.colors.red)
% callbackcalc([],[],hDV)
% end

% decide what sort of hydrology to run
probabilistic1VsDeterministic2= decideDeterministicOrProbHydrology(hDV);


tSlider=get(hDV.hdsldr(3),'val'); % time
switch probabilistic1VsDeterministic2
    case 1 % run proababilistic hydrology with real monthly well data
        %setup well location and rate information from either variable or constant
        %rate structures if running,
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use  %if variable rate injection data loaded
            allWellsDatenumBarrelsPerDay=hDV.data.realWellData.datenumBarrelsPerDay(:,1); % get data for this variable injection rate well
            Xwell=hDV.data.realWellData.XEasting; Ywell=hDV.data.realWellData.YNorthing; %location
        else
            allWellsDatenumBarrelsPerDay=hDV.data.inject.datenumBarrelsPerDay(:,1); % get data for this constant rate well
            Xwell=hDV.data.inject.vals(:,1); Ywell=hDV.data.inject.vals(:,2); %location
        end
        
    case 2 % don't run probabilistic hydro if not is selected, or box x-ed out
        
        if ~isfield(hDV.plotdata.pint,'ppf') || hDV.plotdata.hydprob.yearDataShown ~= tSlider
            calcengine(hDV,'HYDROLOGY')
        end
        %         otherwise make deterministic answer
        % set blue curves to vertical at delta p value
        nfaults=hDV.data.fault.vals(1);
        max_X_data=1;
        for jj=1:nfaults % cycle over faults
            set(hDV.plotdata.dplinesprob(jj),'xdata',[hDV.plotdata.pint.ppf(jj);hDV.plotdata.pint.ppf(jj)],'ydata', [0;1]) % make blue line vertical at deterministic value
            max_X_data=max([hDV.plotdata.pint.ppf(jj),max_X_data]); % find maximum X axis limit
        end
        set(hDV.plotdata.dplinesprob(nfaults+1:hDV.data.NFAULTSMAX),'xdata',[hDV.plotdata.pint.ppf(jj);hDV.plotdata.pint.ppf(jj)],'ydata', [0;1])  % clear old faults off plot
        
        % date on blue vertical axis
        title(hDV.plotdata.hydprob.ax2,{'Pressure on Fault';['Jan 1, ',num2str(tSlider)]}) ;
        hDV.plotdata.hydprob.yearDataShown=tSlider;
        
        set(hDV.plotdata.hydprob.htxtMaxDP,'string',num2str(ceil(max_X_data(1).*1.1)))
        callmaxdpprob(hDV.plotdata.hydprob.htxtMaxDP,[],hDV.plotdata.hydprob.ax2,hDV); % change X axis limit callback function
        set(hDV.plotdata.hydprob.hbutMC,'backgroundcolor',hDV.colors.green)
        return
end

xFault=hDV.data.fault.xf; yFault=hDV.data.fault.yf;
wells = 1:1:hDV.data.nwells; %indices of wells to use
g = hDV.data.adv.vals(6)  ; % constant [m/s^2]

% parameters that can vary
% get storativity and transmissivity
aqthick = hDV.data.reservoir.vals(1)*0.3048 ; %ft to meters
porosityFraction=hDV.data.reservoir.vals(2)/100;
perm_mD= hDV.data.reservoir.vals(3);
rho = hDV.data.adv.vals(5)  ;% [kg/m^3]
dynamicVisc=hDV.data.adv.vals(7) ;%[Pa.s]
FluidCompressibility= hDV.data.adv.vals(8) ;%[Pa^-1]
RockCompressibility= hDV.data.adv.vals(9) ;% [Pa^-1]


indatacell={Xwell,Ywell,wells,xFault,yFault,g,tSlider,aqthick,porosityFraction,perm_mD,rho,...
    dynamicVisc,FluidCompressibility,RockCompressibility,{allWellsDatenumBarrelsPerDay}};
% vary perm and porosity
% sigcell ={0,0,0,0,0,0,0,0,.05,100,0,0,0,0,0} ;
%     nruns =1000 ;

%    probHydrology.sigvals(3)=100;
sigs=hDV.data.probHydrology.sigvals;
sigcell ={0,0,0,0,0,0,0,sigs(1)*0.3048,sigs(2)/100,sigs(3),sigs(4),sigs(5),sigs(6),sigs(7),0}       ;

nruns=hDV.data.probHydrology.sigvals(8);

%monte carlo engine run
% result is pore pressures, 1 column per fault
[outs,hDV.data.probHydrology.distributionsinData] =ppOnFault_MC(indatacell,sigcell,nruns);
hDV.data.probHydrology.outsAllFaults=outs;

max_X_data=1; % X axis max psi
for jj=1:size(outs,2)
    [f,x]=ecdf(outs(:,jj)  );
    probHydrology.fltPresCDFX{jj}=x;
    probHydrology.fltPresCDFY{jj}=1-f; % probability of exceeding that pressure= 1-CDF of pressure
    
    if nargin==0
        figure(1001)
        hold on
        plot(x,1-f)
    else
        set(hDV.plotdata.dplinesprob(jj),'xdata',probHydrology.fltPresCDFX{jj},'ydata', probHydrology.fltPresCDFY{jj})
        max_X_data=max([max(probHydrology.fltPresCDFX{jj}),max_X_data]); % find maximum X axis limit
    end
end

if nargin==0
    figure(1001)
    grid on
    set(figure(1001),'color','w')
    xlabel('pore pressure on fault, PSI')
    ylabel({'probability that that pressure is exceeded';'on that fault'})
else % running in program, not debugging/coding
    set(hDV.plotdata.dplinesprob(hDV.data.fault.vals(1)+1:hDV.data.NFAULTSMAX),'xdata',NaN,'ydata',NaN)
    set(hDV.plotdata.hydprob.hbutMC,'backgroundcolor',hDV.colors.green)
    
    % x axis limit to see blue curves
    if max_X_data==1; max_X_data=str2double(get(hDV.plotdata.pprob.htxtMaxDP,'string'))   ;end % if no data, see CDF curves
    set(hDV.plotdata.hydprob.htxtMaxDP,'string',num2str(ceil(max_X_data(1).*1.1)))
    callmaxdpprob(hDV.plotdata.hydprob.htxtMaxDP,[],hDV.plotdata.hydprob.ax2,hDV); % change X axis limit callback function
    title(hDV.plotdata.hydprob.ax2,{'Probability of Pressure ';['Exceedance on Fault, Jan 1, ',num2str(tSlider)]}) ;
    hDV.plotdata.hydprob.yearDataShown=tSlider;
end


end


%

