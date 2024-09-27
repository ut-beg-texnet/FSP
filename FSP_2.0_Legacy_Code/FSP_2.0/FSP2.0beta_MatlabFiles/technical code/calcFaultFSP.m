%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate FSP
% Given fault number, calculate Fault Slip Potential
% By Rall Walsh, Stanford
% takes probabilistic hydrology, or deterministic

function [FaultFSPs]=calcFaultFSP(hDV,faultstoRun)
if nargin==0 % for debugging
    hDV=evalin('base', 'hSV'); % get hDV to test?
    faultstoRun=1:hDV.data.fault.vals(1);
elseif nargin==1
    faultstoRun=1:hDV.data.fault.vals(1);
    
end


if all(get(hDV.plotdata.pprob.hbutMC,'backgroundcolor')==hDV.colors.red)
    msgWindow1=msgbox('need to run probabilistic if you want to see FSPs', 'calcFaultFSP warning','warn');
    centerFigure(hDV.hfig,msgWindow1);
    return
end
FaultFSPs=zeros(length(faultstoRun),1); % preallocate

for j=faultstoRun
    presCDF = get(hDV.plotdata.flinesprob(j),'xdata');
    probCDF = get(hDV.plotdata.flinesprob(j),'ydata');
    [~,ix]=unique(presCDF); presCDF = presCDF(ix) ; probCDF = probCDF(ix) ;
    
    if length(presCDF)==1 %this is the deterministic geomechanics case
        if hDV.plotdata.pint.ppf(j)<presCDF
            %             hDV.plotdata.pint.fsp(j)=0;
            FaultFSPs(j)=0;
        else
            %             hDV.plotdata.pint.fsp(j)=1;
            FaultFSPs(j)=1;
        end
    else % probabilistic geomech case
        %         FaultFSPs(j,1) = interp1(presCDF,probCDF,hDV.plotdata.pint.ppf(j),'nearest','extrap');         % was for deterministic hydrology
        thisFaultHydroPressures=get(hDV.plotdata.dplinesprob(j),'xdata');
        thisFaultProbsToAverage=zeros(length(thisFaultHydroPressures),1);
        for jj45=1:length(thisFaultHydroPressures)
            thisFaultProbsToAverage(jj45,1)=interp1(presCDF,probCDF,thisFaultHydroPressures(jj45),'nearest','extrap'); % slice Geomec CDF at each blue curve datapoint
            
%             hDV.plotdata.pint.ppf
        end
        FaultFSPs(j,1)=mean(thisFaultProbsToAverage); % take mean of all FSPs of each hydrologic model
    end
end


% hDV.plotdata.pint.fsp=FaultFSPs;


end










