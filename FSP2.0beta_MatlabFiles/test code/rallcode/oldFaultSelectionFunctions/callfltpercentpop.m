%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
function callfltpercentpop(src,~,hDV)
val = get(src,'Value');
if val==1
    set(hDV.plotdata.flinesprob(:), 'visible', 'on');
    set(hDV.plotdata.dplinesprob(:),'visible','on'); % prob. hydro
    set(hDV.plotdata.greyCDFBackground(:),'visible','on'); % prob. hydro
    title(hDV.plotdata.pprob.ax3,'Choose a fault to see sensitivity analysis','fontsize',hDV.ftsz) ;
    title(hDV.plotdata.pflot.ax3,'Mohrs Circles for Fault All','fontsize',12);
    set(hDV.plotdata.pprob.barlow,'ydata',NaN*get(hDV.plotdata.pprob.barlow,'ydata'));
    set(hDV.plotdata.pprob.barhigh,'ydata',NaN*get(hDV.plotdata.pprob.barhigh,'ydata'));
    set(get(hDV.plotdata.pprob.barhigh,'BaseLine'),'BaseValue',0);
    
else
    fon=val-1;
    if 0 % hide curves or grey out?
    set(hDV.plotdata.flinesprob(:), 'visible', 'off');
    set(hDV.plotdata.flinesprob(fon), 'visible', 'on');
    else
    set(hDV.plotdata.flinesprob(:), 'visible', 'off');
    set(hDV.plotdata.flinesprob(fon), 'visible', 'on');  
    end
    title(hDV.plotdata.pprob.ax3,['Sensitivity Analysis for Fault #',num2str(val-1)],'fontsize',hDV.ftsz) ;
    title(hDV.plotdata.pflot.ax3,['Mohrs Circles for Fault ',num2str(val-1)],'fontsize',12);
    set(hDV.plotdata.dplinesprob(:),'visible','off'); % prob. hydro
    set(hDV.plotdata.greyCDFBackground(:),'visible','off'); % prob. hydro
    set(hDV.plotdata.dplinesprob(fon),'visible','on'); % prob. hydro
    set(hDV.plotdata.greyCDFBackground(fon),'visible','on'); % prob. hydro
    
    cl = get(hDV.plotdata.flinesprob(fon),'color') ;
    
    
    names = hDV.plotdata.pprob.names; 
    low_vals = hDV.plotdata.results.barlowdata{fon};
    high_vals = hDV.plotdata.results.barhighdata{fon};
    base_value = hDV.plotdata.results.barnom{fon};
    
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
hDV.plotdata.curfault = val-1 ;
end
