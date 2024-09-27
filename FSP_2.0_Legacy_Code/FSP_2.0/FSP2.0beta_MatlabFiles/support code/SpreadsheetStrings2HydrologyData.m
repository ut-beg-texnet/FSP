% by Rall Walsh
% stanford Summer 2016
% read in text data of X,Y, Pressure, year from uitable
% convert to double precision matrix of data
% could be improved by allowing specification: is it total pressure or
% pressure perturbation?

function SpreadsheetStrings2HydrologyData(stringsImportedHydrology,hDV)



% %%%% only for testing
if nargin ==0 % for coding
 load('./test code\rall code\MODFLOWModelData/examplehDVWithHydrologyData.mat','hDV')
 stringsImportedHydrology=hDV.data.reservoir.stringsImportedHydrology;
end
%%%%% end only for testing

% convert all to doubles from strings
numbersImportedHydrology=str2double(stringsImportedHydrology);
hDV.data.reservoir.numbersImportedHydrology=numbersImportedHydrology;
% find unique years/decimal years
yearsRepresentedHydroImport=unique(hDV.data.reservoir.numbersImportedHydrology(:,4));
hDV.data.reservoir.yearsRepresentedHydroImport=sort(yearsRepresentedHydroImport,'ascend');

%  Set date slider in calculate button?
%  set year limits +-.001 in refreshplotdata? or where? 


if nargin ==0 % for coding
    for index4459=yearsRepresentedHydroImport
    figure(index4459)
    thisTimeData=numbersImportedHydrology(numbersImportedHydrology(:,4)==index4459,:);
    scatter(thisTimeData(:,1),thisTimeData(:,2),30,thisTimeData(:,3),'filled')
    end
end


end % end SpreadsheetStrings2WellData function











