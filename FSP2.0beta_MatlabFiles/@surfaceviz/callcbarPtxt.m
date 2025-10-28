function callcbarPtxt(~,~,hDV)

%get the numbers from GUI
minc = str2double(get(hDV.plotdata.pprob.cminPtxt,'string'));
maxc = str2double(get(hDV.plotdata.pprob.cmaxPtxt,'string'));


%check and reset if needed
if isnan(minc) || isnan(maxc) || minc>=maxc
    %defaults if error
    set(hDV.plotdata.pprob.cminPtxt,'string','0');
    set(hDV.plotdata.pprob.cmaxPtxt,'string','5000');
else
    setprobcolor(hDV);
    
    % update probabilistic colorbar
    mincfc = str2double(get(hDV.plotdata.pprob.cminPtxt,'string'));
    maxcfc = str2double(get(hDV.plotdata.pprob.cmaxPtxt,'string'));
    axis(hDV.plotdata.pprob.axCBAR,[mincfc,maxcfc,0,1]);
    set(hDV.plotdata.pprob.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
    ticksAx4=get(hDV.plotdata.pprob.axCBAR,'xtick');
    ticksAx4(1)=[]; ticksAx4(end)=[];  % don't replot axis labels at ends, where they are set
    set(hDV.plotdata.pprob.axCBAR,'xtick',ticksAx4);
    patchXs=linspace(mincfc,maxcfc, hDV.cmapInfo.patchesInColorbar+1);
    for j=1:hDV.cmapInfo.patchesInColorbar
        % value on colorscale
        fr = (j-1)/(hDV.cmapInfo.patchesInColorbar-1) ; fr2 = (fr.*(maxcfc-mincfc))+mincfc ;
        PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
        colorForPatch = getcolor(flipud(hDV.cmapGYR),fr2,mincfc,maxcfc);
        set(hDV.plotdata.pprob.colorBarPatches(j,1),'facecolor',colorForPatch,'xdata',PatchXCoord)
    end
    
end


end


