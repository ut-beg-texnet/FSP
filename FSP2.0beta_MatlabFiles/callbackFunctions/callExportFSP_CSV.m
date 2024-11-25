% By Rall Walsh, Stanford
% starting in V0.94
% Export data
% could be improved by exporting fault end coordinates, not just middle
% save CDF curve data
function callExportFSP_CSV(~,~,hDV)

[fileNameToPrint,path2] = uiputfile('*.csv','Save Fault Pore Pressure Curve data As'); % find file name and location to save

if any(fileNameToPrint~=0) % if file can exist?
    
    path_file=fullfile(path2,fileNameToPrint);
    TimeData=get(hDV.plotdata.pint.FaultFSPCurvesThruTime(1),'xdata');   % get time data
    fid22 = fopen(path_file,'w');% open file, writing permission
    if fid22== -1 % if couldn't get permission
        errorWindow1=errordlg(cat(2,'couldn''t load or create file: ',path_file,' Maybe it is open? Try closing it so it can be edited. If not, check permissions.  '));
        centerFigure(hDV.hfig,errorWindow1);
        return % don't save
    end
    
    fprintf(fid22, 'fault data:,,1 fault per line,, and 1 date per column ,,,Fault Slip Potential [0 to 1] in same column as corresponding date. Date is January 1st at 12:01 AM in corresponding year:\n');
    fclose(fid22);
    dlmwrite(path_file,TimeData,'-append','coffset',7)
    fid223 = fopen(path_file, 'a'); % append
    fprintf(fid223, ['Fault Number, CenterX, CenterY, Strike(deg), Dip(Deg), length(km),Coefficient of Friction (mu),FSP Curve points (unitless fraction),,\n']);% print headers
    
    fclose(fid223); % close file
    faultsToCycleOver=hDV.plotdata.curfault(1);
    if faultsToCycleOver==0;faultsToCycleOver=[1:1:hDV.data.fault.vals(1)];end % cycle over all faults
    for j443Idx= faultsToCycleOver % cycle over faults to print data
        thisCurve=hDV.plotdata.pint.FaultFSPCurvesThruTime(j443Idx) ;  %this fault
        FSPData=get(thisCurve,'ydata'); % fault FSP in each year
        dlmwrite([path2,fileNameToPrint],[j443Idx,hDV.data.fault.xf(j443Idx),hDV.data.fault.yf(j443Idx),hDV.data.fault.thf(j443Idx),hDV.data.fault.dipf(j443Idx),hDV.data.fault.lenf(j443Idx)...
            ,hDV.data.fault.muf(j443Idx) ,FSPData],'-append');  % ,'coffset',2),RGBDataThisFault(1),RGBDataThisFault(2),RGBDataThisFault(3) ,red (RGB 0to1),green (RGB 0to1),blue (RGB 0to1),
    end
    
end


end

