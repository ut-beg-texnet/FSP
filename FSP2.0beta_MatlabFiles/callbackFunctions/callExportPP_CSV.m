% By Rall Walsh
% starting in V0.94
% Export data
% trouble exporting just 1 fault, unlike exportFSP function
% could be improved by exporting fault end coordinates, not just middle
% save CDF curve data
function callExportPP_CSV(~,~,hDV)

[fileNameToPrint,path2] = uiputfile('*.csv','Save Fault Pore Pressure Curve data As');% find where to put csv created

if any(fileNameToPrint~=0) % if file
    
    path_file=fullfile(path2,fileNameToPrint); % make path
    TimeData=get(hDV.plotdata.pint.PressureCurvesThruTime(1),'xdata');   % get time data
    fid22 = fopen(path_file, 'w');% write permission
    if fid22== -1 % if couldn't get permissions
        errorWindow1=errordlg(cat(2,'couldn''t load or create file: ',path_file,' Maybe it is open? Try closing it so it can be edited. If not, check permissions. ')); % throw error
        centerFigure(hDV.hfig,errorWindow1);
        return
    end
    
    fprintf(fid22, 'fault data:,,1 fault per line,, and 1 date per column ,,,PorePressure (PSI in same column as corresponding date). Date is January 1st at 12:01 AM in corresponding year\n');
    fclose(fid22);% close
    dlmwrite(path_file,TimeData,'-append','coffset',7)% write data
    fid223 = fopen(path_file, 'a'); % append
    fprintf(fid223, ['Fault Number, CenterX, CenterY, Strike(deg), Dip(Deg), length(km),Coefficient of Friction (mu),PorePressure Curve points (in PSI),,\n']); % headers
    fclose(fid223);
    
    faultsToCycleOver=hDV.plotdata.curfault(1);
    if faultsToCycleOver==0;faultsToCycleOver=[1:1:hDV.data.fault.vals(1)];end % cycle over all faults
    for j443Idx= faultsToCycleOver % cycle over faults to print data
        thisCurve=hDV.plotdata.pint.PressureCurvesThruTime(j443Idx) ;        %pick fault curve
        PSIData=get(thisCurve,'ydata');% this fault PSI data
        %                RGBDataThisFault=get(thisCurve,'color'); % when was saving color too
        dlmwrite([path2,fileNameToPrint],[j443Idx,hDV.data.fault.xf(j443Idx),hDV.data.fault.yf(j443Idx),hDV.data.fault.thf(j443Idx),hDV.data.fault.dipf(j443Idx),hDV.data.fault.lenf(j443Idx)...
            ,hDV.data.fault.muf(j443Idx) ,PSIData],'-append');  % ,'coffset',2),RGBDataThisFault(1),RGBDataThisFault(2),RGBDataThisFault(3) ,red (RGB 0to1),green (RGB 0to1),blue (RGB 0to1),
    end% cycle over faults to print data
end % if file
end % function
