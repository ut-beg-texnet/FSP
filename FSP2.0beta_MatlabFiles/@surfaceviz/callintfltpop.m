%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
function callintfltpop(src,~,hDV)
val = get(src,'Value') ;

if hDV.plotdata.printFunctionName
    disp(['running callintfltpop'])
end
ShowAllIntegratedCurves(hDV,val) % if we get the calculation into Calcengine, we can move the contents of ShowAllIntegratedCurves into callintfltpop

if val==1
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
    title(hDV.plotdata.pflot.ax3,'Mohrs Circles for Fault All','fontsize',12);
    
    title(hDV.plotdata.pint.ax4,'FSP vs. Time for All Faults') ;
    set(hDV.plotdata.pint.FaultFSPCurvesThruTime(:), 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
    
else
    fon=val-1;
    set(hDV.plotdata.flinesint(:), 'visible', 'off');
    set(hDV.plotdata.flinesinttxt(:), 'visible', 'off');
    set(hDV.plotdata.pint.FaultFSPCurvesThruTime(:), 'visible', 'off');
    set(hDV.plotdata.flinesint(fon), 'visible', 'on');
    set(hDV.plotdata.flinesinttxt(fon), 'visible', 'on');
    set(hDV.plotdata.pint.FaultFSPCurvesThruTime(fon), 'visible', hDV.plotdata.pint.FSPThruTimeVisibility1); % this makes it so that we can turn the plot off in setupplotpanels
    %     numberOfTimes=100;
    %     tlims= get(hDV.plotdata.pint.ax2,'xlim'); t=linspace(tlims(1) , tlims(2),numberOfTimes);
    tlims= get(hDV.plotdata.pint.ax2,'xlim'); t= round(tlims(1)):1: tlims(2);
    
    %compute pressure as a function of time
    wells=1:1:hDV.data.nwells;
    x=hDV.data.fault.xf(fon);
    y=hDV.data.fault.yf(fon);
    p=pfieldcalc(hDV,x,y,wells,t);
    
    set(hDV.plotdata.pint.fltpressure,'xdata',t,'ydata',p);
    set(hDV.plotdata.pint.fltpressurelim,'xdata',tlims,'ydata',hDV.plotdata.results.outs.ppfail(fon)*[1 1]) ;
    
    c1 = getcolor(hDV.cmapGYR,0.01,hDV.plotdata.minint,hDV.plotdata.maxint);
    c2 = getcolor(hDV.cmapGYR,0.2,hDV.plotdata.minint,hDV.plotdata.maxint);
    c3 = getcolor(hDV.cmapGYR,0.5,hDV.plotdata.minint,hDV.plotdata.maxint);
    %set(hDV.plotdata.pint.backar,'xdata',tlims,'ydata',0*tlims+hDV.plotdata.results.nomhigh{fon},'basevalue',hDV.plotdata.results.nomlow{fon});
    set(hDV.plotdata.pint.fltpressurelimUp,'xdata',tlims,'ydata',hDV.plotdata.results.nomhigh{fon}*[1 1],'displayname',['50% P(slip) at ',num2str(hDV.plotdata.results.nomhigh{fon},6),' PSI'],'color',c3) ;
    set(hDV.plotdata.pint.fltpressurelimDown,'xdata',tlims,'ydata',hDV.plotdata.results.nomlow{fon}*[1 1],'displayname',['20% P(slip) at ',num2str(hDV.plotdata.results.nomlow{fon},6),' PSI'],'color',c2) ;
    set(hDV.plotdata.pint.fltpressurelimBottom,'xdata',tlims,'ydata',hDV.plotdata.results.nomBottom{fon}*[1 1],'displayname',['1% P(slip) at ',num2str(hDV.plotdata.results.nomBottom{fon},6),' PSI'],'color',c1) ;
    
    %                  hDV.plotdata.pint.probHydroMaxMinPp = plot(hDV.plotdata.pint.ax2,NaN,NaN,'b+') ;
    %             hDV.plotdata.pint.probHydroMaxMinPpLabel=text([NaN],[NaN],{'Highest Possible Pp';'Lowest possible Pp'},'parent',hDV.plotdata.pffot.ax2,'color','b','HorizontalAlignment','right',...
    %                 'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','bold','VerticalAlignment','middle','clipping','off');
    
    AllPressuresPlottedYear=get(hDV.plotdata.dplinesprob(fon),'xdata');
    if length(AllPressuresPlottedYear)>2; % deterministic if blue line is only 2 points in probhydro
        set(hDV.plotdata.pint.probHydroMaxMinPp,'xdata',[1;1].*hDV.plotdata.hydprob.yearDataShown,'ydata',[max(AllPressuresPlottedYear);min(AllPressuresPlottedYear)],'visible','on') ;
        set(hDV.plotdata.pint.probHydroMaxMinPpLabel(1),'position',[hDV.plotdata.hydprob.yearDataShown,max(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','bottom','HorizontalAlignment','right')  % label highest and lowest pressure on fault
        set(hDV.plotdata.pint.probHydroMaxMinPpLabel(2),'position',[hDV.plotdata.hydprob.yearDataShown,min(AllPressuresPlottedYear),1],'visible','on','VerticalAlignment','top','HorizontalAlignment','left')
    else
        set(hDV.plotdata.pint.probHydroMaxMinPp,'visible','off') ;
        set(hDV.plotdata.pint.probHydroMaxMinPpLabel,'visible','off') ;
    end
    pmax = max([hDV.plotdata.results.outs.ppfail(fon) ; max(p);max(AllPressuresPlottedYear);1]) ; % max p for ylimit
    set(hDV.plotdata.pint.ax2,'ylim',[0 1.1*pmax]) ;
    
    %edit title
    title(hDV.plotdata.pint.ax2,['Pressure Change at Fault Midpoint for Fault #',num2str(fon)]) ;
    title(hDV.plotdata.pflot.ax3,['Mohrs Circles for Fault ',num2str(fon)],'fontsize',12);
    title(hDV.plotdata.pint.ax4,['FSP vs. Time for Fault #',num2str(fon)]) ;
    %reset time bars to correct height
    set(hDV.plotdata.pint.timebar1,'ydata',get(hDV.plotdata.pint.ax2,'ylim'));
    
    
end
hDV.plotdata.curfault = val-1 ;
end