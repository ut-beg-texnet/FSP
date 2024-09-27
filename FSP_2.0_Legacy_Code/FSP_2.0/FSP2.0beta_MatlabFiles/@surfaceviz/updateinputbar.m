function updateinputbar(hDV)
% next to tornado chart, input variability percentage

if hDV.data.stress.aphi.use==11 % add aphi in tornado chart, remove horizontal stresses
    v1=[hDV.data.sigvals(1) ; hDV.data.sigvals(5:8) ; hDV.data.sigvals(10) ;hDV.data.stress.aphi.sigvals(1)];
    v2=[hDV.data.stress.vals(1) ; hDV.data.stress.vals(6) ; 180 ; 90 ; 180 ; mean(hDV.data.fault.muf);hDV.data.stress.aphi.vals(1)] ;    % order: sV, pp , strike , dip , stress direction , fault friction, aphi
elseif hDV.data.stress.aphi.use==12 % aphi and shmin, no shmax
    v1=[hDV.data.sigvals(1:2) ; hDV.data.sigvals(5:8) ; hDV.data.sigvals(10) ;hDV.data.stress.aphi.sigvals(1)];
    v2=[hDV.data.stress.vals(1:2) ; hDV.data.stress.vals(6) ; 180 ; 90 ; 180 ; mean(hDV.data.fault.muf);hDV.data.stress.aphi.vals(1)] ;    % order: sV, shmin, pp , strike , dip , stress direction , fault friction, aphi
elseif hDV.data.stress.aphi.use==0
    v1=[hDV.data.sigvals(1:3) ; hDV.data.sigvals(5:8) ; hDV.data.sigvals(10) ];
    v2=[hDV.data.stress.vals(1:3) ; hDV.data.stress.vals(6) ; 180 ; 90 ; 180 ; mean(hDV.data.fault.muf)] ;    % order: stress state (3 components), pp , strike , dip , stress direction , fault friction
end


low_vals =-100*v1./v2 ;
high_vals =100*v1./v2;
base_value = 0;
names = hDV.plotdata.pprob.names;

set(hDV.plotdata.pprob2.barlow,'ydata',low_vals,'facecolor',.5*[1 1 1]);
set(hDV.plotdata.pprob2.barhigh,'ydata',high_vals,'facecolor',.5*[1 1 1]);
set(get(hDV.plotdata.pprob2.barhigh,'BaseLine'),'BaseValue',base_value);
set(hDV.plotdata.pprob.ax4,'yticklabel',names);
set(hDV.plotdata.pprob.ax4,'Ytick',1:length(names),'YTickLabel',1:length(names));
set(hDV.plotdata.pprob.ax4,'yticklabel',names);

end
