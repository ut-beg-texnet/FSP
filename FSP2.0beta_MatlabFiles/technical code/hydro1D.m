%  Hydro1D
% by Rall Walsh, Stanford
% make monte-Carlo-Able fuinction that calculates radial
% flow pressure on each fault.
% modeled after Mohrs3D by Darren

function fld=hydro1D(inDataCell)

if nargin==0 % check calculation, test by calling with no input,
    getHDVDataToTest=0; % decide which way to test/debug 1 with GUI or 0 plot 
    if getHDVDataToTest==1 % get data from hDV?
        hDV=evalin('base', 'hSV'); % get hDV to test? 
        %setup well location and rate information from either variable or constant
        %rate structures
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use  %variable rate
            allWellsDatenumBarrelsPerDay=hDV.data.realWellData.datenumBarrelsPerDay(:,1); % get data for this variable injection rate well
            Xwell=hDV.data.realWellData.XEasting; Ywell=hDV.data.realWellData.YNorthing; %location
        else % make constant rate injection into format of variable rate data.
            allWellsDatenumBarrelsPerDay=hDV.data.inject.datenumBarrelsPerDay(:,1); % get data for this constant rate well
            Xwell=hDV.data.inject.vals(:,1); Ywell=hDV.data.inject.vals(:,2); %location
        end
        xFault=hDV.data.fault.xf; yFault=hDV.data.fault.yf;
        tSlider = get(hDV.hdsldr(1),'value') ; % time
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
        
    else % generic inputs. wells in a line
        Xwell=[0.4;10.3];
        Ywell=[0;0];
        % constant well rate data
        allWellsDatenumBarrelsPerDay{1,1}=[   737061.000694444  ,     0;
            749479.000694444  ,       23000];
        allWellsDatenumBarrelsPerDay{2,1}=[   737061.000694444  ,     0;
            749479.000694444  ,       12500];
        xFault=[0:.5:20]';yFault=zeros(size(xFault));
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
        
    end
    % makeIndataCell
inDataCell={Xwell,Ywell,wells,xFault,yFault,g,tSlider,aqthick,porosityFraction,perm_mD,rho,...
dynamicVisc,FluidCompressibility,RockCompressibility,{allWellsDatenumBarrelsPerDay}};
end % end inputs for checking calculation


Xwell=inDataCell{1};
Ywell=inDataCell{2};
wells=inDataCell{3};
xFault=inDataCell{4};
yFault=inDataCell{5};
g=inDataCell{6};
tSlider=inDataCell{7};
aqthick=inDataCell{8};
porosityFraction=inDataCell{9};
perm_mD=inDataCell{10};
rho=inDataCell{11};
dynamicVisc=inDataCell{12};
FluidCompressibility=inDataCell{13};
RockCompressibility=inDataCell{14};
allWellsDatenumBarrelsPerDayCell=inDataCell{15}; % pass as a cell so monte carlo ignores it. 
allWellsDatenumBarrelsPerDay=allWellsDatenumBarrelsPerDayCell{1};

% do calculation
[S , T] = calcST([0,0,0,0,rho,g,dynamicVisc,FluidCompressibility,RockCompressibility]' , aqthick ,porosityFraction,perm_mD) ;
fld = 0.*xFault; % setup data matrix
for index44=1:1:length(wells) %loop over each well to get the field
    
    k=wells(index44);
    xwell = Xwell(k);
    ywell = Ywell(k);
    
    thisWellDatenumBarrelsPerDay=allWellsDatenumBarrelsPerDay{k,1}; % get data for this well
    difflocs=xFault+1i*yFault - (xwell+1i*ywell) ; R=abs(difflocs) ;    %grid distances to well
    rmeters=R*1e3 ; %meters
    
    if  any(year(thisWellDatenumBarrelsPerDay(:,1))<= max(tSlider) ) %only positive time
        pfrontResult = pfront(rmeters,tSlider,thisWellDatenumBarrelsPerDay,S,T,rho,g); % call pfront function
        fld=fld+pfrontResult ; %bars
    end
end % end cycling wells

% psi to bars for consistency, stochastic piece
%smag = hDV.data.reservoir.vals(4)/14.5 ;
%fld = fld+smag*randn(size(fld)) ;

fld(isinf(fld)) = 0 ; %remove infinity values, typically occurs when well centers exactly match a grid location
fld(isnan(fld)) = 0 ; %remove nan values, typically occurs when time at calc is less that the well start times
fld=fld*14.5; %convert to psi
% end

if nargin==0 % code verification
    if getHDVDataToTest==1 % get data from hDV?
        % verify by plotting xs on pressure through time curves in integrated
        % panel
        deletablePlotHandle = plot(tSlider,fld,'parent',hDV.plotdata.pint.ax2,...
            'linewidth',3,'HandleVisibility','off','color','r','marker','x');
    else % generic plot
        figure;
        set(gcf,'color','w')
        plot(xFault,fld,'bx')
        hold on
        plotLimits=get(gca,'ylim');
        plot([Xwell,Xwell]',[plotLimits;plotLimits]','k:')
        xlabel('x coordinate 2D flow of 2 wells in a line')
        ylabel('pressure (PSI) at faults along the line')
    end
end
end


