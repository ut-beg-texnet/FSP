% By Rall Walsh, Stanford
% starting in V0.94
% Export data

function callExportFSPHydrology_CSV(~,~,hDV)

[fileNameToPrint,path2] = uiputfile('*.csv','Save Pore Pressure data As'); % find file name and location to save

if any(fileNameToPrint~=0) % if file can exist?
    
    path_file=fullfile(path2,fileNameToPrint);
    
    % data to export:
    
    hp = hDV.plotdata.pflot ; % hydrology contour
    
    dataToExport=[reshape(hp.Xgrid,[],1),reshape(hp.Ygrid,[],1),reshape(hp.Zgrid,[],1)];
    
    fid22 = fopen(path_file,'w');% open file, writing permission
    if fid22== -1 % if couldn't get permission
        errorWindow1=errordlg(cat(2,'couldn''t load or create file: ',path_file,' Maybe it is open? Try closing it so it can be edited. If not, check permissions.  '));
        centerFigure(hDV.hfig,errorWindow1);
        return % don't save
    end
    ts = get(hDV.hdsldr(1),'value')  ; %years (assume all sliders are identical)
    fprintf(fid22, ['pressure data:,,,,1 coordinate per line,,,,Additional pore pressure. Date is January 1st at 12:01 AM in ',num2str(ts),' \n']);
    fclose(fid22);
%     dlmwrite(path_file,TimeData,'-append','coffset',7)
    fid223 = fopen(path_file, 'a'); % append
    fprintf(fid223, ['X_Easting_km, Y_Northing_km, additionalPressure_PSI,\n']);% print headers
    
    fclose(fid223); % close file
    
    size(dataToExport,1)
    for j443Idx= 1:size(dataToExport,1) % cycle over points to print data
        thisPoint=dataToExport(j443Idx,:);   %this fault
        dlmwrite(path_file,thisPoint,'-append');  %
    end
    
end


end

