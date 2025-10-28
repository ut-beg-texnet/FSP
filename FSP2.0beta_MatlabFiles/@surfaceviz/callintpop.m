%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
%%% for "Choose Plot Labels" in the integrated tab dropdown menu callback


function callintpop(src,~,hDV)
nfaults = hDV.data.fault.vals(1) ; xf = hDV.data.fault.xf ; yf = hDV.data.fault.yf ;

if hDV.plotdata.printFunctionName
    disp(['running callintpop '])
end
% update fault plot numbers
switch get(src,'value')
    case 2
        dv = 1:1:nfaults ; strf = '%10.0f';
    case 3
        dv = hDV.plotdata.pint.ppf ; strf = '%10.0f';
    case 4
        if all(isnan(hDV.plotdata.pint.fsp))
            dv = NaN*(1:1:nfaults); strf = '%10.2f';
        else
            dv = hDV.plotdata.pint.fsp ; strf = '%10.2f';
        end
end

zs = max(hDV.plotdata.pint.Zgrid(:)) ;


if get(src,'value')>1
    plsp=' ';
    if hDV.plotdata.curfault(1)==0 % if all faults selected
        faults=1:nfaults;
    else % one or more faults selected
        faults= hDV.plotdata.curfault;
    end
    for m=faults
        set(hDV.plotdata.flinesinttxt(m), 'string', [plsp,num2str(dv(m),strf)],'position',[xf(m) yf(m) zs],'visible','on');
    end
    
    set(hDV.plotdata.flinesinttxt(nfaults+1:end),'position', [NaN NaN zs]) ;
else
    set(hDV.plotdata.flinesinttxt(:),'string','') ;
    set(hDV.plotdata.flinesinttxt(:),'visible','off') ;
end
callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV,'Integrated')
end