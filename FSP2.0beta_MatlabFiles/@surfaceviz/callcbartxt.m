function callcbartxt(~,~,hDV)

%get the numbers from GUI
mincfc = str2double(get(hDV.plotdata.pffot.cmintxt,'string'));
maxcfc = str2double(get(hDV.plotdata.pffot.cmaxtxt,'string'));


%check and reset if needed
if isnan(mincfc) || isnan(maxcfc) || mincfc>=maxcfc
    set(hDV.plotdata.pffot.cmintxt,'string',num2str(hDV.plotdata.mincfc,'%5.2f'));
    set(hDV.plotdata.pffot.cmaxtxt,'string',num2str(hDV.plotdata.maxcfc));
else
    %update colors 
    hDV.plotdata.mincfc=mincfc;
    hDV.plotdata.maxcfc=maxcfc;
    cv=hDV.plotdata.results.outs.ppfail;
    
    %edit colors for new range
    nfaults = hDV.data.fault.vals(1) ;
    for j=1:1:nfaults
        cl = getcolor(flipud(hDV.cmapGYR),cv(j),hDV.plotdata.mincfc,hDV.plotdata.maxcfc) ;
        set(hDV.plotdata.flinesgeo(j), 'color', cl);
        set(hDV.plotdata.snet(j),'color',cl) ;
        set(hDV.plotdata.snetpoles(j),'markerfacecolor',cl,'markeredgecolor',cl) ;
        set(hDV.plotdata.pffot.mflt(j),'markerfacecolor',cl) ;
        
    end
    
    %for composite plot
    caxis(hDV.plotdata.pffot.ax3,[hDV.plotdata.mincfc,hDV.plotdata.maxcfc]) ;
    
    
    % update  colorbar
    mincfc = str2double(get(hDV.plotdata.pffot.cmintxt,'string'));
    maxcfc = str2double(get(hDV.plotdata.pffot.cmaxtxt,'string'));
    axis(hDV.plotdata.pffot.axCBAR,[mincfc,maxcfc,0,1]);
    set(hDV.plotdata.pffot.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
    ticksAx4=get(hDV.plotdata.pffot.axCBAR,'xtick');
    ticksAx4(1)=[]; ticksAx4(end)=[];  % don't replot axis labels at ends, where they are set
    set(hDV.plotdata.pffot.axCBAR,'xtick',ticksAx4);
    patchXs=linspace(mincfc,maxcfc, hDV.cmapInfo.patchesInColorbar+1);
    for j=1:hDV.cmapInfo.patchesInColorbar
        % value on colorscale
        fr = (j-1)/(hDV.cmapInfo.patchesInColorbar-1) ; fr2 = (fr.*(maxcfc-mincfc))+mincfc ;
        PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
        colorForPatch = getcolor(flipud(hDV.cmapGYR),fr2,mincfc,maxcfc);
        set(hDV.plotdata.pffot.colorBarPatches(j,1),'facecolor',colorForPatch,'xdata',PatchXCoord);
    end
    
    
    
    
end

        %get from current fault value
        refreshFaultSelectorList(hDV,'GEOMECHANICS') % multi fault selector
        callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV)
        
end


