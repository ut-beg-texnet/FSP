function callintcbartxt(~,~,hDV) % integrated colorbar text limit values

%get the numbers from GUI
minint = str2double(get(hDV.plotdata.pint.cmintxt,'string'));
maxint = str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));


%check and reset if needed
if isnan(minint) || isnan(maxint) || minint>=maxint || minint<0 || maxint>1
    set(hDV.plotdata.pint.cmintxt,'string',num2str(hDV.plotdata.minint,'%5.2f'));
    set(hDV.plotdata.pint.cmaxtxt,'string',num2str(hDV.plotdata.maxint));
else
    %update colors
    hDV.plotdata.minint=minint;
    hDV.plotdata.maxint=maxint;
    
    cv=hDV.plotdata.pint.fsp;
    %edit colors for new range
    nfaults = hDV.data.fault.vals(1) ;
    for j=1:1:nfaults
        if ~isnan(cv(j))
            cl = getcolor(hDV.cmapGYR,cv(j),hDV.plotdata.minint,hDV.plotdata.maxint) ;
            set(hDV.plotdata.flinesint(j), 'color', cl);
        else
            set(hDV.plotdata.flinesint(j), 'color', .8*[1 1 1]);
        end
    end
end

mincfc = str2double(get(hDV.plotdata.pint.cmintxt,'string'));
maxcfc = str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));
axis(hDV.plotdata.pint.axCBAR,[mincfc,maxcfc,0,1]);
set(hDV.plotdata.pint.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
ticksAx4=get(hDV.plotdata.pint.axCBAR,'xtick');
ticksAx4(1)=[]; ticksAx4(end)=[];  % don't replot axis labels at ends, where they are set
set(hDV.plotdata.pint.axCBAR,'xtick',ticksAx4);
patchXs=linspace(mincfc,maxcfc, hDV.cmapInfo.patchesInColorbar+1);
for j=1:hDV.cmapInfo.patchesInColorbar
    % value on colorscale
    fr = (j-1)/(hDV.cmapInfo.patchesInColorbar-1) ; fr2 = (fr.*(maxcfc-mincfc))+mincfc ;
    PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
    colorForPatch = getcolor(hDV.cmapGYR,fr2,mincfc,maxcfc);
    set(hDV.plotdata.pint.colorBarPatches(j,1),'facecolor',colorForPatch,'xdata',PatchXCoord)
end
%       

        refreshFaultSelectorList(hDV,'INTEGRATED')
        callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV,'Integrated')

        
end




