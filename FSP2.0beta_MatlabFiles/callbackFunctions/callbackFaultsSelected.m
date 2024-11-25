

% by Rall Walsh
%Stanford

% if one or multiple or all faults selected, run this callback
%  hDV.plotdata.ListboxFaultSelector
%
% callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV,currentTabName)


function callbackFaultsSelected(src,~,hDV,currentTabName)

curfault=get(src,'Val')-1; % currently selected fault numbers as row vector


% get numbers from hDV
% menuHandle=hDV.plotdata.ListboxFaultSelector; % handle of uicontrol object
% nfaults=hDV.data.fault.vals(1); % number of faults
if nargin ==3
    currentTabName=hDV.currtab.name;
end
if any(curfault==0) % if all is selected, even if others are selected too, do only All
    curfault=0;
end
hDV.plotdata.curfault=curfault;


if hDV.plotdata.printFunctionName
    disp(['running callbackFaultsSelected for ',currentTabName])
end

switch currentTabName
    case {'INTEGRATED'}
        % used to be the function  callintfltpop(src,~,hDV)
        
        % ShowAllIntegratedCurves(hDV,val) % if we get the calculation into Calcengine, we can move the contents of ShowAllIntegratedCurves into callintfltpop
        
        if curfault(1)==0 % then all faults
            set(hDV.plotdata.flinesint(:), 'visible', 'on');
            set(hDV.plotdata.flinesinttxt(:), 'visible', 'on');
            set(hDV.plotdata.pint.fltpressure,'xdata',0,'ydata',0);
            set(hDV.plotdata.pint.fltpressurelim,'xdata',0,'ydata',0);
            %set(hDV.plotdata.pint.backar ,'xdata',0,'ydata',0);
            
            set(hDV.plotdata.pint.fltpressurelimUp,'xdata',0,'ydata',0,'displayname','50% FSP','color',[.7 .7 .7]);
            set(hDV.plotdata.pint.fltpressurelimDown,'xdata',0,'ydata',0,'displayname','20% FSP','color',[.7 .7 .7]);
            set(hDV.plotdata.pint.fltpressurelimBottom,'xdata',0,'ydata',0,'displayname','1% FSP','color',[.7 .7 .7]);
            
            set(hDV.plotdata.pint.probHydroMaxMinPp,'visible','off') ;
            set(hDV.plotdata.pint.probHydroMaxMinPpLabel,'visible','off') ;
            
            title(hDV.plotdata.pint.ax2,'Pressure vs. Time for All Faults') ;
            
            title(hDV.plotdata.pint.ax4,'FSP vs. Time for All Faults') ;
            set(hDV.plotdata.pint.FaultFSPCurvesThruTime(:), 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
            set(hDV.plotdata.pint.FSPtextLabelTime2, 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
            set(hDV.plotdata.pint.PressureCurvesThruTime(:),'visible','on')
            
            title(hDV.plotdata.pint.ax4,'All Faults, FSP Through Time','fontsize',hDV.ftsz) ;
            title(hDV.plotdata.pint.ax2,'Select Fault to Plot Pressures','fontsize',hDV.ftsz);
            
        else % one or more but not all faults selected
            set(hDV.plotdata.flinesint(:), 'visible', 'off');
            set(hDV.plotdata.flinesinttxt(:), 'visible', 'off');
            set(hDV.plotdata.pint.FaultFSPCurvesThruTime(:), 'visible', 'off');
            set(hDV.plotdata.pint.PressureCurvesThruTime(:),'visible','off')
            set(hDV.plotdata.pint.FSPtextLabelTime2, 'visible', 'off');
            
            
            
            set(hDV.plotdata.flinesint(curfault), 'visible', 'on');
            set(hDV.plotdata.flinesinttxt(curfault), 'visible', 'on');
            set(hDV.plotdata.pint.FaultFSPCurvesThruTime(curfault), 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
            set(hDV.plotdata.pint.FSPtextLabelTime2(:,curfault), 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
            set(hDV.plotdata.pint.PressureCurvesThruTime(curfault),'visible','on')
            
            if length(curfault)~=1
                str='Selected Faults';
            else
                str=['Fault #',num2str(curfault)];
            end
            title(hDV.plotdata.pint.ax4,['FSP through time for ',str],'fontsize',hDV.ftsz-1) ;
            title(hDV.plotdata.pint.ax2,['Pressure through time on ',str],'fontsize',hDV.ftsz-1);
            
            % set maximum pressure
            
            
            %             % pressure through time y limit and pressure bounds
            %             AllPressuresPlottedYear=get(hDV.plotdata.dplinesprob(curfault),'xdata');
            %             % this got blue curve data
            %             if length(AllPressuresPlottedYear)>2; % deterministic if blue line is only 2 points in probhydro
            %                 set(hDV.plotdata.pint.probHydroMaxMinPp,'xdata',[1;1].*hDV.plotdata.hydprob.yearDataShown,'ydata',[max(AllPressuresPlottedYear);min(AllPressuresPlottedYear)],'visible','on') ;
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel(1),'position',[hDV.plotdata.hydprob.yearDataShown,max(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','bottom','HorizontalAlignment','right')  % label highest and lowest pressure on fault
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel(2),'position',[hDV.plotdata.hydprob.yearDataShown,min(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','top','HorizontalAlignment','left')
            %             else
            %                 set(hDV.plotdata.pint.probHydroMaxMinPp,'visible','off') ;
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel,'visible','off') ;
            %             end
            %             pmax = max([hDV.plotdata.results.outs.ppfail(fon) ; max(p);max(AllPressuresPlottedYear);1]) ; % max p for ylimit
            %             set(hDV.plotdata.pint.ax2,'ylim',[0 1.1*pmax]) ;
            
            
            % fault map labels
            callintpop(hDV.plotdata.pint.pop(1),[],hDV) % refresh fault number/FSP/pore pressure on fault labels
            
            %             %     numberOfTimes=100;
            %             %     tlims= get(hDV.plotdata.pint.ax2,'xlim'); t=linspace(tlims(1) , tlims(2),numberOfTimes);
            %             tlims= get(hDV.plotdata.pint.ax2,'xlim'); t= round(tlims(1)):1: tlims(2);
            %
            %             %compute pressure as a function of time
            %             wells=1:1:hDV.data.nwells;
            %             x=hDV.data.fault.xf(curfault);
            %             y=hDV.data.fault.yf(curfault);
            %             p=pfieldcalc(hDV,x,y,wells,t);
            %
            %             set(hDV.plotdata.pint.fltpressure,'xdata',t,'ydata',p);
            %             set(hDV.plotdata.pint.fltpressurelim,'xdata',tlims,'ydata',hDV.plotdata.results.outs.ppfail(curfault(1))*[1 1]) ;
            %
            %             c1 = getcolor(hDV.cmapGYR,0.01,hDV.plotdata.minint,hDV.plotdata.maxint);
            %             c2 = getcolor(hDV.cmapGYR,0.2,hDV.plotdata.minint,hDV.plotdata.maxint);
            %             c3 = getcolor(hDV.cmapGYR,0.5,hDV.plotdata.minint,hDV.plotdata.maxint);
            %             %set(hDV.plotdata.pint.backar,'xdata',tlims,'ydata',0*tlims+hDV.plotdata.results.nomhigh{fon},'basevalue',hDV.plotdata.results.nomlow{fon});
            %             set(hDV.plotdata.pint.fltpressurelimUp,'xdata',tlims,'ydata',hDV.plotdata.results.nomhigh{fon}*[1 1],'displayname',['50% P(slip) at ',num2str(hDV.plotdata.results.nomhigh{fon},6),' PSI'],'color',c3) ;
            %             set(hDV.plotdata.pint.fltpressurelimDown,'xdata',tlims,'ydata',hDV.plotdata.results.nomlow{fon}*[1 1],'displayname',['20% P(slip) at ',num2str(hDV.plotdata.results.nomlow{fon},6),' PSI'],'color',c2) ;
            %             set(hDV.plotdata.pint.fltpressurelimBottom,'xdata',tlims,'ydata',hDV.plotdata.results.nomBottom{fon}*[1 1],'displayname',['1% P(slip) at ',num2str(hDV.plotdata.results.nomBottom{fon},6),' PSI'],'color',c1) ;
            %
            %             %                  hDV.plotdata.pint.probHydroMaxMinPp = plot(hDV.plotdata.pint.ax2,NaN,NaN,'b+') ;
            %             %             hDV.plotdata.pint.probHydroMaxMinPpLabel=text([NaN],[NaN],{'Highest Possible Pp';'Lowest possible Pp'},'parent',hDV.plotdata.pffot.ax2,'color','b','HorizontalAlignment','right',...
            %             %                 'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','bold','VerticalAlignment','middle','clipping','off');
            %
            %             AllPressuresPlottedYear=get(hDV.plotdata.dplinesprob(fon),'xdata');
            %             if length(AllPressuresPlottedYear)>2; % deterministic if blue line is only 2 points in probhydro
            %                 set(hDV.plotdata.pint.probHydroMaxMinPp,'xdata',[1;1].*hDV.plotdata.hydprob.yearDataShown,'ydata',[max(AllPressuresPlottedYear);min(AllPressuresPlottedYear)],'visible','on') ;
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel(1),'position',[hDV.plotdata.hydprob.yearDataShown,max(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','bottom','HorizontalAlignment','right')  % label highest and lowest pressure on fault
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel(2),'position',[hDV.plotdata.hydprob.yearDataShown,min(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','top','HorizontalAlignment','left')
            %             else
            %                 set(hDV.plotdata.pint.probHydroMaxMinPp,'visible','off') ;
            %                 set(hDV.plotdata.pint.probHydroMaxMinPpLabel,'visible','off') ;
            %             end
            %             pmax = max([hDV.plotdata.results.outs.ppfail(fon) ; max(p);max(AllPressuresPlottedYear);1]) ; % max p for ylimit
            %             set(hDV.plotdata.pint.ax2,'ylim',[0 1.1*pmax]) ;
            %
            %             %edit title
            %             title(hDV.plotdata.pint.ax2,['Pressure Change at Fault Midpoint for Fault #',num2str(fon)]) ;
            %             title(hDV.plotdata.pflot.ax3,['Mohrs Circles for Fault ',num2str(fon)],'fontsize',12);
            %             title(hDV.plotdata.pint.ax4,['FSP vs. Time for Fault #',num2str(fon)]) ;
            %             %reset time bars to correct height
            %             set(hDV.plotdata.pint.timebar1,'ydata',get(hDV.plotdata.pint.ax2,'ylim'));
        end
        
        
    case {'GEOMECHANICS'}
        
        
        % this used to be in the callfltpop function
        if curfault(1)==0 % all faults
            set(hDV.plotdata.flinesgeo(:), 'visible', 'on'); % mapped faults
            set(hDV.plotdata.snet(:), 'visible', 'on'); % stereonet lines
            set(hDV.plotdata.snetpoles(:), 'visible', 'on'); % stereonet poles
            set(hDV.plotdata.pffot.mflt(:), 'visible', 'on'); % mohr
        else % some but not all faults selected
            set(hDV.plotdata.flinesgeo(:), 'visible', 'off'); % mapped faults
            set(hDV.plotdata.pffot.mflt(:), 'visible', 'off'); % mohr
            set(hDV.plotdata.flinesgeo(curfault), 'visible', 'on'); % mapped faults
            set(hDV.plotdata.pffot.mflt(curfault), 'visible', 'on'); % mohr
        end
        
        if curfault(1)==0;curfault=':';end % if all faults selected
        switch get(hDV.plotdata.pffot.popsnet,'value') % dropdown which stereonet selection: composite, normals, lines
            case 1 % fault normals
                set(hDV.plotdata.snet(:), 'visible', 'off');
                set(hDV.plotdata.snetpoles(:), 'visible', 'off');
                set(hDV.plotdata.snetpoles(curfault), 'visible', 'on');
            case 2 % lines
                set(hDV.plotdata.snet(:), 'visible', 'off');
                set(hDV.plotdata.snetpoles(:), 'visible', 'off');
                set(hDV.plotdata.snet(curfault), 'visible', 'on');
            case 3 % composite
                set(hDV.plotdata.snet(:), 'visible', 'off');
                set(hDV.plotdata.snetpoles(:), 'visible', 'off');
                set(hDV.plotdata.snetpoles(curfault), 'visible', 'on');
        end
        callgeopop(hDV.plotdata.pffot.pop(1),[],hDV);
        
        
        
    case {'PROB. GEOMECH','PROB. HYDRO'}
        %  this used to be the callfltpercentpop function
        if curfault(1)==0
            set(hDV.plotdata.flinesprob(:), 'visible', 'on');
            set(hDV.plotdata.dplinesprob(:),'visible','on'); % prob. hydro
            set(hDV.plotdata.greyCDFBackground(:),'visible','on'); % prob. hydro
            title(hDV.plotdata.pprob.ax3,'Choose a fault to see sensitivity analysis','fontsize',hDV.ftsz) ;
            
            set(hDV.plotdata.pprob.barlow,'ydata',NaN*get(hDV.plotdata.pprob.barlow,'ydata'));
            set(hDV.plotdata.pprob.barhigh,'ydata',NaN*get(hDV.plotdata.pprob.barhigh,'ydata'));
            set(get(hDV.plotdata.pprob.barhigh,'BaseLine'),'BaseValue',0);
            
        else % one or more faults, but not all selected
            if 0 % hide curves or grey out?
                set(hDV.plotdata.flinesprob(:), 'visible', 'off');
                set(hDV.plotdata.flinesprob(curfault), 'visible', 'on');
            else
                set(hDV.plotdata.flinesprob(:), 'visible', 'off');
                set(hDV.plotdata.flinesprob(curfault), 'visible', 'on');
            end
            
            title(hDV.plotdata.pprob.ax3,['Sensitivity Analysis for Fault #',num2str(curfault(1))],'fontsize',hDV.ftsz) ;
            
            set(hDV.plotdata.dplinesprob(:),'visible','off'); % prob. hydro
            set(hDV.plotdata.greyCDFBackground(:),'visible','off'); % prob. hydro
            set(hDV.plotdata.dplinesprob(curfault),'visible','on'); % prob. hydro
            set(hDV.plotdata.greyCDFBackground(curfault),'visible','on'); % prob. hydro
            
            % refresh tornado, if multiple faults selected, pick first
            cl = get(hDV.plotdata.flinesprob(curfault(1)),'color') ; % get first fault selected color
            names = hDV.plotdata.pprob.names;
            low_vals = hDV.plotdata.results.barlowdata{curfault(1)};
            high_vals = hDV.plotdata.results.barhighdata{curfault(1)};
            base_value = hDV.plotdata.results.barnom{curfault(1)};
            
            M=[low_vals ; high_vals ; base_value+0*low_vals];
            [~ , ind]=sort(abs(max(M)-min(M)),'ascend');
            
            high_vals=high_vals(ind);
            low_vals=low_vals(ind);
            names=names(ind);
            
            set(hDV.plotdata.pprob.barlow,'ydata',low_vals,'facecolor',cl);
            set(hDV.plotdata.pprob.barhigh,'ydata',high_vals,'facecolor',cl);
            set(get(hDV.plotdata.pprob.barhigh,'BaseLine'),'BaseValue',base_value);
            set(hDV.plotdata.pprob.ax3,'yticklabel',names);
            set(hDV.plotdata.pprob.ax3,'Ytick',1:length(names),'YTickLabel',1:length(names));
            set(hDV.plotdata.pprob.ax3,'yticklabel',names);
        end
        
        
    case 'HYDROLOGY'
        
        hp = hDV.plotdata.pflot ; % hydrology contour plot
        % find if all, 1 or some faults selected
        if curfault(1)==0
            title(hDV.plotdata.pflot.ax3,'Mohr Circles for All Faults ','fontsize',12);
            faults=1:hDV.data.fault.vals(1);
        elseif length(curfault)==1% % if only 1 fault selected
            title(hDV.plotdata.pflot.ax3,['Mohr Circles for Fault #',num2str(curfault(1))],'fontsize',12);
            faults=hDV.plotdata.curfault;
        else % multiple faults selected
            faults=hDV.plotdata.curfault;
            title(hDV.plotdata.pflot.ax3,['Mohr Circles for ',num2str(length(curfault)),' selected Faults'],'fontsize',12);
        end
        
        % turn all faults off
        set(hp.mflt(:),'visible','off') ;
        set(hp.mc1(:),'visible','off');
        set(hp.mc2(:),'visible','off');
        set(hp.mc3(:),'visible','off');
        
        for j=faults  % show selected faults
            set(hp.mflt(j),'visible','on') ;
            set(hp.mc1(j),'visible','on');
            set(hp.mc2(j),'visible','on');
            set(hp.mc3(j),'visible','on');
        end
        
    otherwise  % case 'MODEL INPUTS'
        % could make fault bigger or label or something
end % end switch tab name




end % end function


