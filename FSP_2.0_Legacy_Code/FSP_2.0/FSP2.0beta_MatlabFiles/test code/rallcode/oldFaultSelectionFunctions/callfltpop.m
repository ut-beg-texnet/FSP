%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
function callfltpop(src,~,hDV)
val = get(src,'Value') ;
if val==1
    set(hDV.plotdata.flinesgeo(:), 'visible', 'on');
    set(hDV.plotdata.snet(:), 'visible', 'on');
    set(hDV.plotdata.snetpoles(:), 'visible', 'on');
    set(hDV.plotdata.pffot.mflt(:), 'visible', 'on');
    title(hDV.plotdata.pflot.ax3,'Mohrs Circles for Fault All','fontsize',12);
else
    fon=val-1;
    set(hDV.plotdata.flinesgeo(:), 'visible', 'off');
    set(hDV.plotdata.snet(:), 'visible', 'off');
    set(hDV.plotdata.snetpoles(:), 'visible', 'off');
    set(hDV.plotdata.pffot.mflt(:), 'visible', 'off');
    
    set(hDV.plotdata.flinesgeo(fon), 'visible', 'on');
    set(hDV.plotdata.snet(fon), 'visible', 'on');
    set(hDV.plotdata.snetpoles(fon), 'visible', 'on');
    set(hDV.plotdata.pffot.mflt(fon), 'visible', 'on');
    title(hDV.plotdata.pflot.ax3,['Mohrs Circles for Fault ',num2str(fon)],'fontsize',12);

    
end
hDV.plotdata.curfault = val -1;
callgeopop(hDV.plotdata.pffot.pop(1),[],hDV);
end