% make histograms of input data on a given fault
% by Rall Walsh
% Stanford 2016
% could be improved by bringing in actual bounds, not just looking at data
% bounds
% this is a callback function to an export button
% modified by Rall Outside Stanford
%
% Modified by Suvrat Lele March 2018 for FSP 2.0; changes related to A-Phi
% Subsurface Mechanics and Induced Seismicity, Drilling & Subsurface
% ExxonMobil URC

function makeInputHistogramsProbabilistic(~,~,hDV)

if nargin==0
    hDV=evalin('base', 'hSV'); % get hDV to test?
end

curfault=hDV.plotdata.curfault(1);
if curfault==0 % if a fault is not selected, throw this
    msgWindow1=msgbox({cat(2,'Select a fault to see it''s data!')}, 'Displaying Inputs','warn');
    centerFigure(hDV.hfig,msgWindow1);
    return
end

% make figure
hfig = hDV.hfig;
hpos = get(hfig,'position') ;
fignew = figure('color',[1 1 1],'units','pixels','visible','off','Name',['Inputs For Fault Number ',num2str(curfault)],'NumberTitle','off');
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
inputtedData=hDV.data.distributionsinData{curfault,1};
depthFt=hDV.data.stress.vals(5);
% column headers
%     indatacell = {Sig0   ,ignoreParameter4 , p0 , strike , dip , SHdir, dp , mu ,
%     biot , nu} ; % towo more optional are aphi and mu

variableNames={'S_{Vertical}';'S_{hmin}';'S_{Hmax}';'Natural Pore Pressure';['Strike of fault #',num2str(curfault)];['Dip of fault #',num2str(curfault)];'S_{Hmax} Azimuth';'Mu'};
unitLabels={'[PSI/ft]';'[PSI/ft]';'[PSI/ft]';'[PSI/ft]';'[Degrees]';'[Degrees]';'[Degrees]';'Coefficient of Friction'};
% Sh,Sv,Sh,InitialPp,Strike,Dip,SHdir,mu
threeStresses=cat(1,inputtedData{:,1})./depthFt;

% j=ecdfhist(get(hDV.plotdata.flinesprob(curfault),'Ydata')',get(hDV.plotdata.flinesprob(curfault),'Xdata')');
dataMatrix=[threeStresses,[inputtedData{:,3}]'./depthFt,[inputtedData{:,4}]',[inputtedData{:,5}]',[inputtedData{:,6}]',[inputtedData{:,8}]'   ];



% hDV.data.distributions.deterministicVals % this is a cell of strings
% txt2 = {'Vertical Stress Grad [',' psi/ft]' ;
%     'Min Horiz. Grad [',' psi/ft]' ;
%     'Max Horiz. Grad [',' psi/ft]' ;
%     'IGNORE PARAMETER4: ' ,'';
%     'Initial PP Grad [',' psi/ft]' ;
%     'Strike Angles [',' degrees]' ;
%     'Dip Angles [',' degrees]' ;
%     'Max Horiz. Stress Dir [',' degrees]' ;
%     'delta pore pressure','' ;
%     'Friction Coeff Mu [',']' ;
%     'Biot Parameter' ,'';
%     'Poisson''s Ratio [' ,']'};
%     'APhi' is 13th and is optional
noshow = [4 9 11 12] ; %what not to show

%find analytical distribution points
DistributionSigvals=hDV.data.sigvals;
DistributionSigvals(noshow)=[];


thf=hDV.data.fault.thf(curfault)  ;   %fault strikes
dips=hDV.data.fault.dipf(curfault) ; %fault dips
SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
mufs=hDV.data.fault.muf(curfault) ;

deterministicVals=[hDV.data.stress.vals(1,1);hDV.data.stress.vals(2,1);hDV.data.stress.vals(3,1);hDV.data.stress.vals(6);thf;dips;SHdir;mufs];

if hDV.data.stress.aphi.use
    dataMatrix=[dataMatrix,[inputtedData{:,11}]'];
    deterministicVals=[deterministicVals;hDV.data.stress.aphi.vals(1)];
    DistributionSigvals=[DistributionSigvals;hDV.data.stress.aphi.sigvals(1)];
    variableNames=[variableNames;{'aPhi parameter'}];
    unitLabels=[unitLabels;{' '}];
end

loValue=deterministicVals-DistributionSigvals;
hiValue=deterministicVals+DistributionSigvals;


% add ppf result to matrix: 
loValue=[loValue;NaN];
hiValue=[hiValue;NaN];
ppfs = get(hDV.plotdata.flinesprob(curfault),'xdata')';
dataMatrix=[dataMatrix,ppfs(2:end)];
variableNames=[variableNames;{['result: pore pressure to slip for fault number ',num2str(curfault)]}];
unitLabels=[unitLabels;{'[PSI]'}];



numPlots=sum(~all(~diff(dataMatrix)));
numrows=2;
numcols=ceil(numPlots./numrows);
indexesToPlot=1:size(dataMatrix,2);
indexesToPlot=indexesToPlot(~all(~diff(dataMatrix)));
nbins=25;
VerticalCoord=size(dataMatrix,1)/nbins; % geomechanics requires 1000 bootstraps, hydrology you can pick number of bootstraps
for plotIndex=1:numPlots
    
    thisAxisIndex=indexesToPlot(plotIndex);
    thisPlotCol=rem( plotIndex,numcols);
    if thisPlotCol==0; thisPlotCol=numcols; end
    thisPlotRow=ceil(thisAxisIndex./numcols);
    
    subplotIndex=plotIndex;
    if mod(plotIndex,2)==1 && plotIndex==max(numPlots);subplotIndex=[plotIndex,plotIndex+1];end
    axisToPlot(plotIndex,1)=subplot(numrows,numcols,subplotIndex,'parent',fignew);
    hold(axisToPlot(plotIndex,1),'on')
   
    
    hist(axisToPlot(plotIndex,1),dataMatrix(:,thisAxisIndex),nbins,'edgecolor','none','facecolor',hDV.colors.blue) % ,'EdgeColor','w','width',1,'FaceColor',[0,0,0.2]
    
%        if plotIndex<=7
    plotHandle22(plotIndex)= plot(axisToPlot(plotIndex,1),[loValue(thisAxisIndex);loValue(thisAxisIndex);hiValue(thisAxisIndex);hiValue(thisAxisIndex)],...
        [0;VerticalCoord;VerticalCoord;0],'r-','linewidth',2.5);
%        end
    
       
    box on
    axis(axisToPlot(plotIndex,1),'tight')
    
    xlims=get(axisToPlot(plotIndex,1),'xlim');
    ylims=get(axisToPlot(plotIndex,1),'ylim');
    ylims=[ylims(1), ylims(2).*1.1];
    set(axisToPlot(plotIndex,1),'ylim',ylims)
    text(xlims(1),ylims(2),num2str(xlims(1)),'verticalalignment','top','horizontalalignment','left')
    text(xlims(2),ylims(2),num2str(xlims(2)),'verticalalignment','top','horizontalalignment','right')
    set(axisToPlot(plotIndex,1),'ytick',[])
    title(axisToPlot(plotIndex,1),variableNames(thisAxisIndex))
    xlabel(unitLabels(thisAxisIndex))
    if thisPlotCol==1
    ylabel(axisToPlot(plotIndex,1),{'{\color[rgb]{1 0 0} Probability Distribution Inputted';['\color[rgb]{',num2str(hDV.colors.blue),'} Number of Realizations}']})
    end
end

end













