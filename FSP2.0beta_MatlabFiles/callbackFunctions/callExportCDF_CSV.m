
% save CDF curve data
% this is called for 
function callExportCDF_CSV(src,~,hDV,curvehandles)

[fileNameToPrint,path2] = uiputfile('*.csv','Save Fault CDF Curve data As');

if any(fileNameToPrint~=0)
    
     path_file=fullfile(path2,fileNameToPrint);
% fileID1=fopen(path_file,'w');

ProbabilityData=get(curvehandles(1),'ydata');   % get probability data
% ProbabilityData=get(hDV.plotdata.flinesprob(1),'ydata');   % get probability data
fid22 = fopen(path_file, 'w');
fprintf(fid22, ['# fault data:,,',get(src,'TooltipString'),',,,,,Probability (PSI in corresponding column below):\n']);
fclose(fid22);
% dlmwrite(path_file,{'Probabilities (Of PSI in this column):','And this'})
dlmwrite(path_file,ProbabilityData,'-append','coffset',7)
fid223 = fopen(path_file, 'a');
fprintf(fid223, ['# Fault Number, CenterX, CenterY, Strike(deg), Dip(Deg), length(km),Coefficient of Friction (mu),CDF Curve points (in PSI), 1 fault per line, and 1 probability per column \n']);

% last is fault number ',num2str(hDV.data.fault.vals(1)),'Probabilities (Of PSI in corresponding column):
fclose(fid223);

% dlmwrite([path2,fileNameToPrint],['Fault CDF Curves in PSI: 1 fault per line, first is fault number 1, last is fault number ',num2str(hDV.data.fault.vals(1))],'-append')
if hDV.plotdata.curfault(1)==0
    faults=1:1:hDV.data.fault.vals(1);
else
    faults=hDV.plotdata.curfault;
end
for j443Idx=faults
%                thisCurve= hDV.plotdata.flinesprob(j443Idx) ;        % =plot(0,0,'linewidth',2,'color','k','parent',hDV.plotdata.pprob.ax2) ;
                thisCurve= curvehandles(j443Idx) ;   
               PSIData=get(thisCurve,'xdata');
%                RGBDataThisFault=get(thisCurve,'color');

               dlmwrite([path2,fileNameToPrint],[j443Idx,hDV.data.fault.xf(j443Idx),hDV.data.fault.yf(j443Idx),hDV.data.fault.thf(j443Idx),hDV.data.fault.dipf(j443Idx),hDV.data.fault.lenf(j443Idx)...
                 ,hDV.data.fault.muf(j443Idx) ,PSIData],'-append');  % ,'coffset',2),RGBDataThisFault(1),RGBDataThisFault(2),RGBDataThisFault(3) ,red (RGB 0to1),green (RGB 0to1),blue (RGB 0to1),
            end



% fclose(fileID1)

end


end












