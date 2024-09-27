% callRunSummaryPlot
% By Rall
% make Fault PP and Fault FSP curves through time

% run for all faults
% calculate deterministic and if necessary probabilistic hydrology through
% time

% calculate deterministic and FSP and probabilistic Hi/lo bounds through
% time

% run for whole duration of model for now.

function callRunSummaryPlot(src,~,hDV)

% % if calculate button red, run that first
% if all(get(hDV.bCalc,'backgroundcolor')==hDV.colors.red)
% callbackcalc([],[],hDV)
% end

startyear=get(hDV.hdsldr(1),'min');
endyear=get(hDV.hdsldr(1),'max'); % this is hardwired in slider bar
nfaults = hDV.data.fault.vals(1); % number of faults
% probabilistic1VsDeterministic2= decideDeterministicOrProbHydrology(hDV);

% find which years to calculate
if hDV.data.reservoir.importHydrology==1
    YearsToRun=hDV.data.reservoir.yearsRepresentedHydroImport;
else
    YearsToRun=startyear:1:endyear;
end

% preallocate
SavedFSPmatrix=zeros(nfaults,length(YearsToRun));
PressuresAllFaultsAndTimes=SavedFSPmatrix;
minFaultPP=SavedFSPmatrix;
maxFaultPP=SavedFSPmatrix;


handleWaitbar1=waitbar(0.01,'Calculating FSP January 1 each year','CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
centerFigure(hDV.hfig,handleWaitbar1)
setappdata(handleWaitbar1,'canceling',0)
runCancelled=0;


for kjj449=1:length(YearsToRun)
    waitbar(kjj449./length(YearsToRun),handleWaitbar1,['Calculating FSP in ',num2str(YearsToRun(kjj449))]);
    pause(.3) % .5 second
    set(hDV.hdsldr(3),'val',YearsToRun(kjj449)) % update waitbar value
    sldrcall(hDV.hdsldr(3),[],hDV) % run slider calculation, this calls all probabilisitic calculations necessary
    if getappdata(handleWaitbar1,'canceling') % if cancel button is hit, stop
        disp('Summary plot generation cancelled')
        runCancelled=1;
        break
    end
    % save FSP values this year
    SavedFSPmatrix(:,kjj449)=hDV.plotdata.pint.fsp;
    PressuresAllFaultsAndTimes(:,kjj449) = hDV.plotdata.pint.ppf; % save deterministic pore pressure value
    for j=1:nfaults
        % get lo and high values of fault pore pressure perturbation from
        % blue curve
        % if deterministic, these will be the same as the det. value
        minFaultPP(j,kjj449) = min(get(hDV.plotdata.dplinesprob(j),'xdata'));
        maxFaultPP(j,kjj449) = min(get(hDV.plotdata.dplinesprob(j),'xdata'));
    end
    
end
delete(handleWaitbar1)  % DELETE the waitbar; don't try to CLOSE it.

% call make plot function
if ~runCancelled
    set(src,'backgroundcolor',hDV.colors.green)
    %     PlotSummaryPlot(SavedFSPmatrix,YearsToRun,hDV,faultsSelected,valsout(3))
end


% refresh background color patches
FSPYlims=get(hDV.plotdata.pint.ax4,'ylim');
patchYs=linspace(FSPYlims(1),FSPYlims(2), length(hDV.plotdata.pint.FSPBackgroundPatches)+1); % Y coordinates of bounds of patches
PatchXCoord= [startyear,endyear,endyear,startyear];% X coordinates of patches
minint = str2double(get(hDV.plotdata.pint.cmintxt,'string'));
maxint = str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));
for j=1:length(hDV.plotdata.pint.FSPBackgroundPatches)% cycle over each background patch
    PatchYCoord = [patchYs(j),patchYs(j),patchYs(j+1),patchYs(j+1)] ;
    cl = getcolor(hDV.cmapGYR,mean([patchYs(j),patchYs(j+1)]),minint,maxint) ; % color at this FSP value, RGB, get color at middle of patch
    set(hDV.plotdata.pint.FSPBackgroundPatches(j),'xdata',PatchXCoord,'ydata',PatchYCoord,...
        'facecolor',cl,'visible',hDV.plotdata.pint.FSPThruTimeVisibility1);
end

% update plots
for jj7=1:1:nfaults   % set fault positions - each fault is a row of PressuresAllFaultsAndTimes
    set(hDV.plotdata.pint.PressureCurvesThruTime(jj7),'xdata',YearsToRun,'ydata',PressuresAllFaultsAndTimes(jj7,:))
    set(hDV.plotdata.pint.FaultFSPCurvesThruTime(jj7),'xdata',YearsToRun,'ydata',SavedFSPmatrix(jj7,:))
    
    
end

if nfaults<hDV.data.NFAULTSMAX         % remove old faults
    set(hDV.plotdata.pint.PressureCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
    set(hDV.plotdata.pint.FaultFSPCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
end

maxValueForYlim=max([max(max(PressuresAllFaultsAndTimes)),1]); % get largest value to find axis limits
set(hDV.plotdata.pint.ax2,'ylim',[0,1.05*maxValueForYlim])    % set Y limit of Axis

%reset green vertical time bar to correct height
set(hDV.plotdata.pint.timebar1,'ydata',get(hDV.plotdata.pint.ax2,'ylim'));


ticks=get(hDV.plotdata.pint.ax4,'xtick');
ticks=ticks(ticks==round(ticks)); % remove non whole years
set(hDV.plotdata.pint.ax4,'xtick',ticks);
% hDV.plotdata.pint.FSPtextLabelTime2=cell(nfaults,1);
numTextLabels=size(hDV.plotdata.pint.FSPtextLabelTime2,1);
set(hDV.plotdata.pint.FSPtextLabelTime2,'position',[NaN,NaN]); % reset all text label positions to off plot
for j=1:nfaults
    FSPLineLabelYcoord=[];
    FSPLineLabelXcoord=[];
    
    % decide which points to label
    for kindx445=1:length(YearsToRun) % cycle over year index
        [uniqueValues,faultsWithUniqueValues,~]=unique(SavedFSPmatrix(:,kindx445),'first'); % find first instance of each FSP value
        if any(faultsWithUniqueValues==j) && uniqueValues(faultsWithUniqueValues==j)>0.0011 %*(max(FSPYlims)-min(FSPYlims)) % then we want a label at this location and year
            % also make a requirement that value is greater than 1%?
            FSPLineLabelYcoord(end+1,1)=uniqueValues(faultsWithUniqueValues==j); % save label coordinates
            FSPLineLabelXcoord(end+1,1)=YearsToRun(kindx445);
        end
        
        
        
        %             % downsample if lots of years, try to stagger too
        %             if length(FSPLineLabelXcoord)>=41 && ~isempty(FSPLineLabelXcoord)
        %                 FSPLineLabelXcoord=FSPLineLabelXcoord(mod(j,5)+1:8:end,1);
        %                 FSPLineLabelYcoord=FSPLineLabelYcoord(mod(j,5)+1:8:end,1);
        %             elseif length(FSPLineLabelXcoord)>=21 && ~isempty(FSPLineLabelXcoord)
        %                 FSPLineLabelXcoord=FSPLineLabelXcoord(mod(j,3)+1:4:end,1);
        %                 FSPLineLabelYcoord=FSPLineLabelYcoord(mod(j,3)+1:4:end,1);
        %             elseif length(FSPLineLabelXcoord)>=6 && ~isempty(FSPLineLabelXcoord)
        %                 FSPLineLabelXcoord=FSPLineLabelXcoord(mod(j,2)+1:2:end,1);
        %                 FSPLineLabelYcoord=FSPLineLabelYcoord(mod(j,2)+1:2:end,1);
        %             end
        if length(FSPLineLabelXcoord)>numTextLabels
            indexes=randperm(length(FSPLineLabelXcoord)); % reshuffle indexes
            indexes=indexes(1:numTextLabels); % take first 10 reshuffled indexes
            indexes=sort(indexes);
            FSPLineLabelYcoord=FSPLineLabelYcoord(indexes);
            FSPLineLabelXcoord=FSPLineLabelXcoord(indexes);
        end
        
        if ~isempty(FSPLineLabelXcoord)
            for  jj55=1:length(FSPLineLabelXcoord) % cycle over text labels
                
                set(hDV.plotdata.pint.FSPtextLabelTime2(jj55,j),'position',[FSPLineLabelXcoord(jj55),FSPLineLabelYcoord(jj55)])
                %         this line below is in setupplotpanels
                %          hDV.plotdata.pint.FSPtextLabelTime2(1:faultTextLabelsMax,j)= text(nan(faultTextLabelsMax,1),nan(faultTextLabelsMax,1),num2str(j),...
                %                    'verticalalignment','bottom','color','b','parent',hDV.plotdata.pint.ax4);
                
            end
        end
    end
end



% % for refreshing faults
% %             if ~all(gret(hDV.plotdata.pprob.summaryPlotButton,'backgroundcolor')==hDV.colors.red)
% %                 set(vertcat(hDV.plotdata.pint.FSPtextLabelTime2{:}),'visible','on')
% %             else
% %                set(vertcat(hDV.plotdata.pint.FSPtextLabelTime2{:}),'visible','off')
% %             end

% if  isfield(hDV.plotdata,'results') && strcmpi(hDV.currtab.name,'integrated') % if all faults selected, and probabilistic has been run, and on integrated panel
%     % if internal hydrology calculation, not imported hydrology calculation
%     if hDV.data.reservoir.importHydrology==0 % if using well data, not imported hydrologic model
%         tlims=[startyear,endyear]; timesOnPlot= round(tlims(1)):1: tlims(2);   % time years of interest
%         PressuresAllFaultsAndTimes=zeros(nfaults,length(timesOnPlot)); % preallocate
%         fspAllFaultsAndTimes=PressuresAllFaultsAndTimes; % preallocate
%         %wait bar
%         hWaitbar = waitbar(0.01,'Calculating ...','Name','FSP through time');
%         centerFigure(hDV.hfig,hWaitbar)
%
%         for fon= 1:1:length(timesOnPlot );
%             if ~ishandle(hWaitbar)
%                 disp('All faults FSP Premature Shutdown') ;
%                 break
%             end
%             wells=1:1:hDV.data.nwells; % wells to use: all
%             x=hDV.data.fault.xf; %  % find all fault coordinates
%             y=hDV.data.fault.yf ;%
%             PressuresAllFaultsAndTimes(:,fon)=pfieldcalc(hDV,x,y,wells,timesOnPlot(fon)); % get pressure on all faults at this time
%             waitbar(fon / length( timesOnPlot ),hWaitbar,['Calculating ',num2str(round(100*fon / length( timesOnPlot ))), '% ...']) ;
%         end
%         if ishandle(hWaitbar) % close waitbar
%             close(hWaitbar) ;
%         end
%
%
%     else %if hDV.data.reservoir.importHydrology==1 % if imported hydrologic model, like from modflow
%         %wait bar
%         hWaitbar = waitbar(0.01,'Calculating ...','Name','FSP through time');
%         centerFigure(hDV.hfig,hWaitbar)
%
%         % overwrite timesOnPlot
%         timesOnPlot=hDV.data.reservoir.yearsRepresentedHydroImport;
%         PressuresAllFaultsAndTimes=zeros(nfaults,length(timesOnPlot)); % preallocate
%         fspAllFaultsAndTimes=PressuresAllFaultsAndTimes; % preallocate
%         tlims=[min(timesOnPlot),max(timesOnPlot)];
%         %         if tlims(1)==tlims(2); tlims(2)=tlims(2)+0.001; end
%
%         for fon= 1:1:length(timesOnPlot );
%             if ~ishandle(hWaitbar)
%                 disp('All faults FSP Premature Shutdown') ;
%                 break
%             end
%             % interpolate pressure on faults
%             ts = timesOnPlot(fon);
%             numbersImportedHydrology=hDV.data.reservoir.numbersImportedHydrology;
%
%             thisYearData=numbersImportedHydrology(numbersImportedHydrology(:,4)==ts,:);
%             thisYearX=thisYearData(:,1);
%             thisYearY=thisYearData(:,2);
%             thisYearPSI=thisYearData(:,3);
%
%             % interpolate imported hydrology onto faults
%             PressuresAllFaultsAndTimes(:,fon) =griddata(thisYearX,thisYearY,thisYearPSI,hDV.data.fault.xf,hDV.data.fault.yf) ;
%
%             waitbar(fon / length( timesOnPlot ),hWaitbar,['Calculating ',num2str(round(100*fon / length( timesOnPlot ))), '% ...']) ;
%
%         end
%
%         if ishandle(hWaitbar) % close waitbar
%             close(hWaitbar) ;
%         end
%
%     end
%
%     probabilistic1VsDeterministic2= decideDeterministicOrProbHydrology(hDV)
%
%



%         if strcmpi(hDV.plotdata.pint.FSPThruTimeVisibility1,'on') % this can be gotten rid of if we keep (or get rid of FSP through time plot
%             PressuresCDFThisFault = get(hDV.plotdata.flinesprob(jj7),'Xdata');
%             probababilitiesCDFThisFault=get(hDV.plotdata.flinesprob(jj7),'Ydata');
%             %         disp(['fault number ',num2str(jj7),' updating pressure'])
%             for k77=1:1:length( timesOnPlot )
%                 thisFaultAndTimePSI=PressuresAllFaultsAndTimes(jj7,k77);
%                 % this line is outdated now that probabilistic hydrology is
%                 % added:
%                 fspAllFaultsAndTimes(jj7,k77)=max(probababilitiesCDFThisFault(PressuresCDFThisFault==max([max(PressuresCDFThisFault(PressuresCDFThisFault<thisFaultAndTimePSI)),min(PressuresCDFThisFault)])));
%             end
%             disp('incorrect FSP calculation in ShowAllIntegratedCurves (outdated by probabilistic hydrology)')
%             set(hDV.plotdata.pint.FaultFSPCurvesThruTime(jj7),'xdata',timesOnPlot,'ydata',fspAllFaultsAndTimes(jj7,:))
%         end

% end

% end






















