% by Rall Walsh
% Stanford August 2016
% plot all faults Pressure, Time curves, colored by trafficlight value
% experimented with lineChangingColor.m
% might be more efficient if it loops over times, not faults, because
% faults are vectorized in pfront, but time isn't. 
% maybe even plot by fault loop, but calculate by time loop. 
% was called from callintfltpop.m
% these lines are in setupplotpanels.m
%             for j=1:1:hDV.data.NFAULTSMAX
%                hDV.plotdata.pint.ColoredPressureCurvesThruTime(j) = surface([NaN,NaN;NaN,NaN],[NaN,NaN;NaN,NaN],[NaN,NaN;NaN,NaN],[NaN,NaN;NaN,NaN],'parent',hDV.plotdata.pint.ax2,...
%                'facecol','no','edgecol','interp','linew',2,'HandleVisibility','off');
%              end


function ColoredIntegratedCurves(hDV,dropDownMenuValue)


tlims= get(hDV.plotdata.pint.ax2,'xlim'); timesOnPlot= ceil(tlims(1)):1: tlims(2);   % time years of interest
nfaults = hDV.data.fault.vals(1); % number of faults
if dropDownMenuValue==1  &&  isfield(hDV.plotdata,'results') && strcmpi(hDV.currtab.name,'integrated')% if all faults selected, and probabilistic has been run,

    %compute pressure as a function of time
    wells=1:1:hDV.data.nwells; % wells to use: all
    maxValueForYlim=1; % reset axis ylim (PSI)
    
%     cMinProbabilisticColorbar=str2double(get(hDV.plotdata.pint.cmintxt,'string'));
%     cMaxProbabilisticColorbar =  str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));
%     if isnan(cMaxProbabilisticColorbar);cMaxProbabilisticColorbar=1;end
%     if isnan(cMinProbabilisticColorbar);cMinProbabilisticColorbar=0;end
    
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
        
        % make background trafficLightColored
%         PressuresCDFThisFault = get(hDV.plotdata.flinesprob(fon),'Xdata');
%         probababilitiesCDFThisFault=get(hDV.plotdata.flinesprob(fon),'Ydata');
        
%         plotColorRGB=zeros(2,length(timesOnPlot),3); % clear/preallocate
%         for jjkk7=1:length(timesOnPlot) % build C data
%             thisPatchPSIToColorWith=PressuresThroughTimeThisFault(jjkk7); %
%             % find highest pressure in this CDF curve below this fault PSI value
%             thisPointProbabilityToColorWith=max(probababilitiesCDFThisFault(PressuresCDFThisFault==max([max(PressuresCDFThisFault(PressuresCDFThisFault<thisPatchPSIToColorWith)),min(PressuresCDFThisFault)])));
%             
%             thisPointColor=getcolor(hDV.cmapGYR,thisPointProbabilityToColorWith,cMinProbabilisticColorbar,cMaxProbabilisticColorbar);
%             plotColorRGB(1:2,jjkk7,1)=thisPointColor(1); % RGB on 3 pages
%             plotColorRGB(1:2,jjkk7,2)=thisPointColor(2);
%             plotColorRGB(1:2,jjkk7,3)=thisPointColor(3);
%         end
        
        % surface objects need matrix inputs, so repeat the line data twice
%         plotTimeData=[timesOnPlot;timesOnPlot]; % time data two rows repeated
%         plotPressureData=[PressuresThroughTimeThisFault;PressuresThroughTimeThisFault]; % pressure data 2 rows repeated
%         plotZdata=zeros(size(plotPressureData)); % z coord = 0
%         plotColorRGB=ones(size(plotPressureData)).*0.5;
        
               set(hDV.plotdata.pint.PressureCurvesThruTime(fon),'xdata',timesOnPlot,'ydata',PressuresThroughTimeThisFault)

           
%         set(hDV.plotdata.pint.ColoredPressureCurvesThruTime(fon),'xdata',plotTimeData,'ydata',plotPressureData,'zdata',plotZdata,'cdata',plotColorRGB); % set line values
        maxValueForYlim=max([max(PressuresThroughTimeThisFault),maxValueForYlim]); % save largest value to find axis limits
        waitbar(fon / nfaults,hWaitbar,['Calculating ',num2str(100*fon / nfaults), '% ...']) ;
        
    end
    % remove old faults
    set(hDV.plotdata.pint.PressureCurvesThruTime(nfaults+1:end),'xdata',NaN,'ydata',NaN)
%     set(hDV.plotdata.pint.ColoredPressureCurvesThruTime(nfaults+1:end),'xdata',[NaN,NaN;NaN,NaN],'ydata',[NaN,NaN;NaN,NaN],'zdata',[NaN,NaN;NaN,NaN],'cdata',[NaN,NaN;NaN,NaN]);
    set( hDV.plotdata.pint.ax2,'ylim',[0,1.1*maxValueForYlim])    % set Y limit of Axis
    if ishandle(hWaitbar) % close waitbar
        close(hWaitbar) ;
    end
    
else % dropdown menu selected a given fault, so hide the all faults curves
       set(hDV.plotdata.pint.PressureCurvesThruTime(:),'xdata',NaN,'ydata',NaN)
%     set(hDV.plotdata.pint.ColoredPressureCurvesThruTime(:),'xdata',[NaN,NaN;NaN,NaN],'ydata',[NaN,NaN;NaN,NaN],'zdata',[NaN,NaN;NaN,NaN],'cdata',[NaN,NaN;NaN,NaN]);
end

if ~isfield(hDV.plotdata,'results') && any([strcmpi(hDV.currtab.name,'integrated'),strcmpi(hDV.currtab.name,'hydrology')])
    edlgBox=errordlg('Run Probabilistic Tab First!','Data Error') ;centerFigure(hDV.hfig,edlgBox);
    
end

end


