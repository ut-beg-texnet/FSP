% make histograms of input data on a given fault
% by Rall Walsh
% Stanford 2016
% could be improved by bringing in actual bounds, not just looking at data
% bounds
% this is a callback function to an export button

 function callMakeInputHistogramsHydro(~,~,hDV)

if nargin==0
    hDV=evalin('base', 'hSV'); % get hDV to test?
end

% if probabilistic hydrology hasn't been run at all, or is red. 
if ~isfield(hDV.data.probHydrology,'distributionsinData') ||  all(get(hDV.plotdata.hydprob.hbutMC,'backgroundcolor')==hDV.colors.red)
    msgWindow1=msgbox({cat(2,'run probabilistic hydrology analysis first!')}, 'Displaying Inputs','warn');
    centerFigure(hDV.hfig,msgWindow1);
    return
end

%find analytical distribution points
DistributionSigvals=hDV.data.probHydrology.sigvals;


% make figure
hfig = hDV.hfig;
hpos = get(hfig,'position') ;
fignew = figure('color',[1 1 1],'units','pixels','visible','off','Name',['Uncertain Hydrologic Inputs: ',num2str(DistributionSigvals(8)),' Realizations'],'NumberTitle','off');
pos = [0 0 .9*hpos(3) .6*hpos(4)];
set(fignew,'position',pos) ;
set(fignew,'MenuBar','none');
set(fignew,'NumberTitle','off');
set(fignew,'DefaultUicontrolUnits','normalized');
set(fignew,'DefaultUicontrolFontsize',14);
% set(fignew,'closerequestfcn',@crf) ;
centerFigure(hfig,fignew);
% imshow('help.jpg','border','tight');
set(fignew,'visible','on') ;

% get data for this fault
inputtedData=hDV.data.probHydrology.distributionsinData;
% depthFt=hDV.data.stress.vals(5);
% column headers

%     indatacell={Xwell,Ywell,wells,xFault,yFault,g,tSlider,aqthick,porosityFraction,perm_mD,rho,...
%         dynamicVisc,FluidCompressibility,RockCompressibility,{allWellsDatenumBarrelsPerDay}};

%   % parameters that can vary
%     % get storativity and transmissivity
%     aqthick = 100*0.3048 ; %100 ft to 30.48 meters
%     porosityFraction=0.1;
%     perm_mD= 200;
%     rho = 1000 ; %fluid density [kg/m^3]
%     dynamicVisc= 0.0008 ;%[Pa.s]
%     FluidCompressibility=  3.6e-10 ;%[Pa^-1]
%     RockCompressibility=  1.08e-09 ;% [Pa^-1]

variableNames={'Injection Formation Thickness';'Porosity';'Permeability';'Fluid Density';'Dynamic Viscosity';'Fluid Compressibility';'Rock Compressibility'};
unitLabels={'[ft]';'[%]';'[mD]';'[kg/m^3]';'[Pa.s]';'[Pa^{-1}]';'[Pa^{-1}]'};
% Sh,Sv,Sh,InitialPp,Strike,Dip,SHdir,mu
% threeStresses=cat(1,inputtedData{:,1})./depthFt;

dataMatrix= [[inputtedData{:,8}]'./0.3048 ,[inputtedData{:,9}]'.*100,[inputtedData{:,10}]',[inputtedData{:,11}]',[inputtedData{:,12}]',[inputtedData{:,13}]',[inputtedData{:,14}]'];   % [threeStresses,[inputtedData{:,3}]'./depthFt,[inputtedData{:,4}]',[inputtedData{:,5}]',[inputtedData{:,6}]',[inputtedData{:,8}]' ];


% find current fault data
if hDV.plotdata.curfault(1) ~=0 
set(fignew,'name',['Hydrologic Uncertainties on Fault Number ',num2str(hDV.plotdata.curfault(1)),' in the Year ',num2str(hDV.plotdata.hydprob.yearDataShown),', (',num2str(DistributionSigvals(8)),' Realizations). '])
variableNames{end+1}=['Result: Pressure on Fault #',num2str(hDV.plotdata.curfault(1)),' in the Year ',num2str(hDV.plotdata.hydprob.yearDataShown)];
unitLabels{end+1}='[PSI]';
dataMatrix(:,end+1)=hDV.data.probHydrology.outsAllFaults(:,hDV.plotdata.curfault(1));
end



% parameters that can vary
% get storativity and transmissivity
aqthick = hDV.data.reservoir.vals(1);%*0.3048 ; %ft to meters
porosityFraction=hDV.data.reservoir.vals(2);% /100;
perm_mD= hDV.data.reservoir.vals(3);
rho = hDV.data.adv.vals(5)  ;% [kg/m^3]
dynamicVisc=hDV.data.adv.vals(7) ;%[Pa.s]
FluidCompressibility= hDV.data.adv.vals(8) ;%[Pa^-1]
RockCompressibility= hDV.data.adv.vals(9) ;% [Pa^-1]

% calculate deterministic, low and high values for red curve inputs
det=[aqthick,porosityFraction,perm_mD,rho,dynamicVisc,FluidCompressibility,RockCompressibility]';
loValue=  det-DistributionSigvals(1:7);
hiValue=  det+DistributionSigvals(1:7);

numPlots=sum(~all(~diff(dataMatrix))); % find number of columns that aren't all equal
numrows=2;
numcols=ceil(numPlots./numrows);
histogramHandles=zeros(numPlots,1);
plotHandle22=zeros(numPlots,1);
nBins=25;
VerticalCoord=DistributionSigvals(8)./nBins; % height of uniform distribution - y coordinate
indexesToPlot=1:size(dataMatrix,2);
  indexesToPlot=indexesToPlot(~all(~diff(dataMatrix)));
for cycles44=1:numPlots
    plotIndex=indexesToPlot(cycles44); % variable number, which might not be same as plot number
    
    thisPlotCol=rem( cycles44,numcols);
    if thisPlotCol==0; thisPlotCol=numcols; end
    thisPlotRow=ceil(cycles44./numcols);
    
    subplotIndex=cycles44;
    if mod(numPlots,2)==1 && plotIndex==8;subplotIndex=[cycles44,cycles44+1];end
    axisToPlot(cycles44,1)=subplot(numrows,numcols,subplotIndex,'parent',fignew);
    hold(axisToPlot(cycles44,1),'on')
    
    
    hist(axisToPlot(cycles44,1),dataMatrix(:,plotIndex),nBins,'EdgeColor','none','facecolor',hDV.colors.blue) % ,'EdgeColor','w','width',1,'FaceColor',[0,0,0.2]
    if plotIndex<=7
    plotHandle22= plot(axisToPlot(cycles44,1),[loValue(plotIndex);loValue(plotIndex);hiValue(plotIndex);hiValue(plotIndex)],...
        [0;VerticalCoord;VerticalCoord;0],'r-','linewidth',2.5);
    end
    box on
      axis(axisToPlot(cycles44,1),'tight')
    
    xlims=get(axisToPlot(cycles44,1),'xlim');
    ylims=get(axisToPlot(cycles44,1),'ylim');
    ylims=[ylims(1), ylims(2).*1.1];
    set(axisToPlot(cycles44,1),'ylim',ylims)
    if plotIndex==8  % if PSI, limit to tenths of a PSI
        text(xlims(1),ylims(2),num2str(xlims(1),'%.1f'),'verticalalignment','top','horizontalalignment','left')
        text(xlims(2),ylims(2),num2str(xlims(2),'%.1f'),'verticalalignment','top','horizontalalignment','right')
    else
        text(xlims(1),ylims(2),num2str(xlims(1)),'verticalalignment','top','horizontalalignment','left')
        text(xlims(2),ylims(2),num2str(xlims(2)),'verticalalignment','top','horizontalalignment','right')
    end
    
    set(axisToPlot(cycles44,1),'ytick',[])
    title(axisToPlot(cycles44,1),variableNames(plotIndex))
    xlabel(unitLabels(plotIndex))
    if thisPlotCol==1
    ylabel(axisToPlot(cycles44,1),{'{\color[rgb]{1 0 0} Probability Distribution Inputted';['\color[rgb]{',num2str(hDV.colors.blue),'} Number of Realizations}']})
    end
end

 end













