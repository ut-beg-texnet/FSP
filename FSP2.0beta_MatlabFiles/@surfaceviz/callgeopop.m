%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
% geomechanics label selection callback
function callgeopop(src,~,hDV)
nfaults = hDV.data.fault.vals(1) ; xf = hDV.data.fault.xf ; yf = hDV.data.fault.yf ;
% update fault plot numbers
switch get(src,'value')
    case 2
        dv = 1:1:nfaults ; strf = '%10.0f';
    case 3
        dv = hDV.data.fault.dipf ; strf = '%10.0f';
    case 4
        dv = hDV.data.fault.thf ; strf = '%10.0f';
    case 5
        dv = hDV.plotdata.results.outs.ppfail ; strf = '%10.0f';
        dv(dv<0) = 0;  % Display negative pressures to slip as zero
    case 6
        dv = hDV.plotdata.results.outs.cff ; strf = '%10.0f';
    case 7
        dv = hDV.plotdata.results.outs.scu ; strf = '%10.2f';
end

if get(src,'value')>1
    plsp=' ';
    if  hDV.plotdata.curfault(1)==0 % all faults selected
        for m=1:1:nfaults
            set(hDV.plotdata.flinesgeotxt(m), 'string', [plsp,num2str(dv(m),strf)],'position',[xf(m) yf(m)],'visible','on');
        end
    else % one or more faults selected
        set(hDV.plotdata.flinesgeotxt(:),'visible','off','FontSize',12);
        for m=hDV.plotdata.curfault
            set(hDV.plotdata.flinesgeotxt(m), 'string', [plsp,num2str(dv(m),strf)],'position',[xf(m) yf(m)],'visible','on');
        end
    end
    set(hDV.plotdata.flinesgeotxt(nfaults+1:end),'position', [NaN NaN]) ;
    uistack(hDV.plotdata.flinesgeotxt(:),'top');
else
    set(hDV.plotdata.flinesgeotxt(:),'visible','off') ;
end
end