function resetprobpanel(hDV)

set(hDV.plotdata.flinesprob(:),'xdata',0,'ydata',0);
set(hDV.plotdata.pprob.hbutMC,'backgroundcolor','red') ; 
hDV.plotdata.pint.ppm=NaN*hDV.plotdata.pint.ppm;

end