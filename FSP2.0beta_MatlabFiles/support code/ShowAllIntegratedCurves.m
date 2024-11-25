% by Rall Walsh
% Stanford August 2016
% plot all faults Pressure, Time curves,
% was called from callintfltpop.m
% these lines are in setupplotpanels.m
%              for j=1:1:hDV.data.NFAULTSMAX
%                hDV.plotdata.pint.PressureCurvesThruTime(j) = plot(0,0,'parent',hDV.plotdata.pint.ax2,...
%                'linewidth',2,'HandleVisibility','off','color','b');
%              end
% it is now quicker because it puts all faults in at once - the hydrology
% is vectorized in faults but loops in time.
% it's likely we want to calculate the PressuresAllFaultsAndTimes matrix somewhere
% else (calcengine) save it in hDV. Then we can just show it here, or even in callintfltpop. (we should talk about this).

function ShowAllIntegratedCurves(hDV,dropDownMenuValue) % version 2 playing with adding FSP through time plot

nfaults = hDV.data.fault.vals(1); % number of faults

if hDV.plotdata.printFunctionName
    disp(['running ShowAllIntegratedCurves'])
end

if dropDownMenuValue==1  &&  isfield(hDV.plotdata,'results') && strcmpi(hDV.currtab.name,'integrated') % if all faults selected, and probabilistic has been run, and on integrated panel
    
    % if internal hydrology calculation, not imported hydrology calculation
    if hDV.data.reservoir.importHydrology==0 % if using well data, not imported hydrologic model
        tlims=[get(hDV.hdsldr(1),'min'),get(hDV.hdsldr(1),'max')]; timesOnPlot= round(tlims(1)):1: tlims(2);   % time years of interest
        PressuresAllFaultsAndTimes=zeros(nfaults,length(timesOnPlot)); % preallocate
        fspAllFaultsAndTimes=PressuresAllFaultsAndTimes; % preallocate
        %wait bar
        hWaitbar = waitbar(0,'Calculating ...','Name','FSP through time');
        centerFigure(hDV.hfig,hWaitbar)
        
        for fon= 1:1:length(timesOnPlot );
            if ~ishandle(hWaitbar)
                msgWindow1=msgbox(['All faults FSP through time Premature Shutdown'], 'premature shutdown','warn');
                centerFigure(hDV.hfig,msgWindow1);
                break
            end
            wells=1:1:hDV.data.nwells; % wells to use: all
            x=hDV.data.fault.xf; %  % find all fault coordinates
            y=hDV.data.fault.yf ;%
            PressuresAllFaultsAndTimes(:,fon)=pfieldcalc(hDV,x,y,wells,timesOnPlot(fon)); % get pressure on all faults at this time
            waitbar(fon / length( timesOnPlot ),hWaitbar,['Calculating ',num2str(round(100*fon / length( timesOnPlot ))), '% ...']) ;
            %         disp('a')
        end
        if ishandle(hWaitbar) % close waitbar
            close(hWaitbar) ;
        end
        
    else %if hDV.data.reservoir.importHydrology==1 % if imported hydrologic model, like from modflow
        %wait bar
        hWaitbar = waitbar(0,'Calculating ...','Name','FSP through time');
        centerFigure(hDV.hfig,hWaitbar)
        
        % overwrite timesOnPlot
        timesOnPlot=hDV.data.reservoir.yearsRepresentedHydroImport;
        PressuresAllFaultsAndTimes=zeros(nfaults,length(timesOnPlot)); % preallocate
        fspAllFaultsAndTimes=PressuresAllFaultsAndTimes; % preallocate
        tlims=[min(timesOnPlot),max(timesOnPlot)];
        %         if tlims(1)==tlims(2); tlims(2)=tlims(2)+0.001; end
        
        for fon= 1:1:length(timesOnPlot );
            if ~ishandle(hWaitbar)
                msgWindow1=msgbox(['All faults FSP Premature Shutdown'], 'premature shutdown','warn');
                centerFigure(hDV.hfig,msgWindow1);
                break
            end
            % interpolate pressure on faults
            ts = timesOnPlot(fon);
            numbersImportedHydrology=hDV.data.reservoir.numbersImportedHydrology;
            
            thisYearData=numbersImportedHydrology(numbersImportedHydrology(:,4)==ts,:);
            thisYearX=thisYearData(:,1);
            thisYearY=thisYearData(:,2);
            thisYearPSI=thisYearData(:,3);
            
            % interpolate imported hydrology onto faults
            PressuresAllFaultsAndTimes(:,fon) =griddata(thisYearX,thisYearY,thisYearPSI,hDV.data.fault.xf,hDV.data.fault.yf) ;
            
            waitbar(fon / length( timesOnPlot ),hWaitbar,['Calculating ',num2str(round(100*fon / length( timesOnPlot ))), '% ...']) ;
            
        end
        
        if ishandle(hWaitbar) % close waitbar
            close(hWaitbar) ;
        end
        
    end
    
    
    
    % update plots
    for jj7=1:1:nfaults   % set fault positions - each fault is a row of PressuresAllFaultsAndTimes
        set(hDV.plotdata.pint.PressureCurvesThruTime(jj7),'xdata',timesOnPlot,'ydata',PressuresAllFaultsAndTimes(jj7,:))
        
        if strcmpi(hDV.plotdata.pint.FSPThruTimeVisibility1,'on') % this can be gotten rid of if we keep (or get rid of FSP through time plot
            PressuresCDFThisFault = get(hDV.plotdata.flinesprob(jj7),'Xdata');
            probababilitiesCDFThisFault=get(hDV.plotdata.flinesprob(jj7),'Ydata');
            %         disp(['fault number ',num2str(jj7),' updating pressure'])
            for k77=1:1:length( timesOnPlot )
                thisFaultAndTimePSI=PressuresAllFaultsAndTimes(jj7,k77);
                % this line is outdated now that probabilistic hydrology is
                % added:
                fspAllFaultsAndTimes(jj7,k77)=max(probababilitiesCDFThisFault(PressuresCDFThisFault==max([max(PressuresCDFThisFault(PressuresCDFThisFault<thisFaultAndTimePSI)),min(PressuresCDFThisFault)])));
            end
            disp('incorrect FSP calculation in ShowAllIntegratedCurves (outdated by probabilistic hydrology)')
            set(hDV.plotdata.pint.FaultFSPCurvesThruTime(jj7),'xdata',timesOnPlot,'ydata',fspAllFaultsAndTimes(jj7,:))
        end
        
    end
    if nfaults<hDV.data.NFAULTSMAX         % remove old faults
        set(hDV.plotdata.pint.PressureCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
        set(hDV.plotdata.pint.FaultFSPCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
    end
    maxValueForYlim=max([max(max(PressuresAllFaultsAndTimes)),1]); % get largest value to find axis limits
    set( hDV.plotdata.pint.ax2,'ylim',[0,1.05*maxValueForYlim])    % set Y limit of Axis
    
    %reset green vertical time bar to correct height
    set(hDV.plotdata.pint.timebar1,'ydata',get(hDV.plotdata.pint.ax2,'ylim'));
    
else % dropdown menu selected a given fault, so hide the all faults curves
    set(hDV.plotdata.pint.PressureCurvesThruTime(:),'xdata',NaN,'ydata',NaN)
end

end


