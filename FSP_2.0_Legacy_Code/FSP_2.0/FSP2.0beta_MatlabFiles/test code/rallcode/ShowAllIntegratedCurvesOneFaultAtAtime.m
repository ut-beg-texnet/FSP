% by Rall Walsh
% Stanford August 2016
% plot all faults Pressure, Time curves,
% might be more efficient if it loops over times, not faults, because
% faults are vectorized in pfront, but time isn't.
% maybe even plot by fault loop, but calculate by time loop.
% was called from callintfltpop.m
% these lines are in setupplotpanels.m
%              for j=1:1:hDV.data.NFAULTSMAX
%                hDV.plotdata.pint.PressureCurvesThruTime(j) = plot(0,0,'parent',hDV.plotdata.pint.ax2,...
%                'linewidth',2,'HandleVisibility','off','color','b');
%              end

function ShowAllIntegratedCurves(hDV,dropDownMenuValue)


tlims= get(hDV.plotdata.pint.ax2,'xlim'); timesOnPlot= ceil(tlims(1)):1: tlims(2);   % time years of interest
nfaults = hDV.data.fault.vals(1); % number of faults
if dropDownMenuValue==1  &&  isfield(hDV.plotdata,'results') && strcmpi(hDV.currtab.name,'integrated')% if all faults selected, and probabilistic has been run, and on integrated panel
    
    %compute pressure as a function of time
    wells=1:1:hDV.data.nwells; % wells to use: all
    maxValueForYlim=1; % reset axis ylim (PSI)
    
    %wait bar
    hWaitbar = waitbar(0,'Calculating ...','Name','FSP through time');
    centerFigure(hDV.hfig,hWaitbar)
    
    for fon=1:1:nfaults
        
        if ~ishandle(hWaitbar)
            disp('All faults FSP Premature Shutdown') ;
            break
        end
        
        x=hDV.data.fault.xf(fon) ; % find this fault coordinates
        y=hDV.data.fault.yf(fon);
        PressuresThroughTimeThisFault=pfieldcalc(hDV,x,y,wells,timesOnPlot); % get pressure on this fault throgh time
        
        set(hDV.plotdata.pint.PressureCurvesThruTime(fon),'xdata',timesOnPlot,'ydata',PressuresThroughTimeThisFault)
        maxValueForYlim=max([max(PressuresThroughTimeThisFault),maxValueForYlim]); % save largest value to find axis limits
        waitbar(fon / nfaults,hWaitbar,['Calculating ',num2str(100*fon / nfaults), '% ...']) ;
        
    end
    % remove old faults
    set(hDV.plotdata.pint.PressureCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
    set( hDV.plotdata.pint.ax2,'ylim',[0,1.1*maxValueForYlim])    % set Y limit of Axis
    if ishandle(hWaitbar) % close waitbar
        close(hWaitbar) ;
    end
else % dropdown menu selected a given fault, so hide the all faults curves
    set(hDV.plotdata.pint.PressureCurvesThruTime(:),'xdata',NaN,'ydata',NaN)
end

end


