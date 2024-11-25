function setstressregtext(hDV)

%handles to set 
h =[hDV.plotdata.pffot.srtxt;hDV.plotdata.inputMap.srtxt];

if hDV.data.stress.aphi.use==0
    sig = hDV.data.stress.vals(1:3); sig=sort(sig,'descend'); %in situ stress
    ixSv = find(sig==hDV.data.stress.vals(1)); %index of vertical stress
    ixSv = ixSv(1); %use only one of the indices if there is more than one (edge cases)

    switch ixSv
        case 1
            set(h(:),'string','Stress Regime: Normal Faulting');
        case 2
            set(h(:),'string','Stress Regime: Strike-Slip Faulting');
        case 3
            set(h(:),'string','Stress Regime: Reverse Faulting');
    end
else
    if hDV.data.stress.aphi.vals(1) <=1
        set(h(:),'string','Stress Regime: Normal Faulting');
    elseif hDV.data.stress.aphi.vals(1) <=2
        set(h(:),'string','Stress Regime: Strike-Slip Faulting');
    elseif hDV.data.stress.aphi.vals(1) <=3
        set(h(:),'string','Stress Regime: Reverse Faulting');
    end
end


end