% By Rall
% make uimenu to export all possible figures
% could improve by adding .eps files?
function makeExportUImenu(hDV)
% hDV=hSV;
mainfig=hDV.hfig;

% tabnames for each axis that can be exported
axisTitles={'fault and well map';'Well Rates';'fault map';'mohr diagram';'stereonet';...
    'CDF curves';'Variability of Inputs';'Fault Sensitivity Analysis';'Hydrology Map';'Single Well Radial Solutions';...
    'Hydrologic Mohr circles';'fault pressure vs probability';'fault FSP map';'Pressure through time';'FSP through time'};

axisHandles=[hDV.plotdata.inputMap.ax1;hDV.plotdata.inputMap.ax2;hDV.plotdata.pffot.ax;hDV.plotdata.pffot.ax2;hDV.plotdata.pffot.ax3;... % last on this line is stereonet
    hDV.plotdata.pprob.ax2;hDV.plotdata.pprob.ax4;hDV.plotdata.pprob.ax3;hDV.plotdata.pflot.ax;hDV.plotdata.pflot.ax2;hDV.plotdata.pflot.ax3;... % last was hydrologic mohr circles
    hDV.plotdata.hydprob.ax2;hDV.plotdata.pint.ax;hDV.plotdata.pint.ax2;hDV.plotdata.pint.ax4];

hDV.uiMenuHandles.ExportFigMenu.Tophandle(1)=uimenu(mainfig,'Label','Export Image');
hDV.uiMenuHandles.ExportFigMenu.tabnames=zeros(length(hDV.tabnames)+1,1);
axislevelExportUimenus=zeros(length(axisTitles),1);
numberOfAxesOnEachTab = [2;3;3;3;1;3];
axisCount=0;

for k=1:length(hDV.tabnames) % cycle over tabs
    
    hDV.uiMenuHandles.ExportFigMenu.tabnames(k,1)=uimenu(hDV.uiMenuHandles.ExportFigMenu.Tophandle,'Label',hDV.tabnames{k});
    
    for jj5=1:numberOfAxesOnEachTab(k) % cycle over axes
        axisCount=axisCount+1; % count number of axes (to find axisTitle and handle)
        uimenu(hDV.uiMenuHandles.ExportFigMenu.tabnames(k,1),'Label',axisTitles{axisCount},'callback',{@exportFigure,hDV,axisCount,axisTitles,axisHandles});
        
    end
    
end

hDV.uiMenuHandles.ExportFigMenu.tabnames(k+1,1)=uimenu(hDV.uiMenuHandles.ExportFigMenu.Tophandle,'Label','Change Export Image Settings','separator','on');
hDV.uiMenuHandles.ExportFigMenu.settingsMenu =   uimenu(hDV.uiMenuHandles.ExportFigMenu.tabnames(k+1,1),'Label','resolution');
DPIOptions=[150,300,450,600,900,1200];
checks={'off','off','on','off','off','off'};
for j338=1:length(DPIOptions)
    hDV.uiMenuHandles.ExportFigMenu.resolutionSelectorMenus(j338,1) = uimenu(hDV.uiMenuHandles.ExportFigMenu.settingsMenu,'Label',[num2str(DPIOptions(j338)),' dots per inch'],'callback',{@exportFigureTypeSettings},...
        'checked',checks{j338});
end

% can select filetype with uimenu now
% hDV.uiMenuHandles.ExportFigMenu.fileTypeMenu =   uimenu(hDV.uiMenuHandles.ExportFigMenu.tabnames(k+1,1),'Label','file type');
% fileTypeOptions={'jpeg';'png';'tiff'};
% checks2={'on';'off';'off'};
% for j339=1:length(fileTypeOptions)
%     hDV.uiMenuHandles.ExportFigMenu.fileSelectorMenus(j339,1)=uimenu(hDV.uiMenuHandles.ExportFigMenu.fileTypeMenu,'Label',fileTypeOptions{j339},'callback',{@exportFigureTypeSettings},...
%         'checked',checks2{j339});
% end

end

% export selected axis
function exportFigure(~,~,hDV,axisCount,axisTitles,axisHandles)

% bring window to find where to save file
[fileNameToPrint,path2] = uiputfile({'*.jpg';'*.png';'*.tif';'*.eps';'*.fig'},['Save ',axisTitles{axisCount},' Image']);
% if no file selected, don't save
if ~any(fileNameToPrint~=0) % if no file selected, return
    return
end

% combine filename with pathname
path_file=fullfile(path2,fileNameToPrint); % path to save
% desired axis to save
handleax=axisHandles(axisCount);
% make dummy figure
if isdeployed % don't show figure at all if not coding right now
    newFigToDelete=figure('visible','off');
else % show briefly
    newFigToDelete=figure;
end

set(newFigToDelete,'PaperPositionMode','auto')
newax=copyobj(handleax,newFigToDelete);% copy axis to figure
set(newax,'outerposition',[0,0,1,1])% fill figure with axis
% 
if all(get(handleax,'color')==hDV.colors.axsBgrnd) % if background is colored grey so you can see yellow faults
    set(newFigToDelete,'inverthardcopy','off')
end
if strcmpi('stereonet',axisTitles{axisCount}) % stereonet colormap fix bug
    set(newFigToDelete,'colormap',flipud(hDV.cmapGYR))
end
rendererString = ['-',lower(get(hDV.hfig,'Renderer'))]; % get renderer type from main figure (edited in advanced tab
% % find file extension from menu
% checkedFileBinary=strcmpi(get(hDV.uiMenuHandles.ExportFigMenu.fileSelectorMenus,'checked'),'on'); % find which menu option selected
% fileExtensions=get(hDV.uiMenuHandles.ExportFigMenu.fileSelectorMenus,'label'); %jpeg, tiff, etc
% extensionToUse=fileExtensions(checkedFileBinary);% find single string selected
switch path_file(end-3:end) % last 4 characters
    case '.jpg'
        extensionToUse='-djpeg';
    case '.tif'
        extensionToUse='-dtiff';
    case '.png'
        extensionToUse='-dpng';
    case '.eps'
        extensionToUse='-depsc';
    case '.fig'
        set(newFigToDelete,'visible','on')
        saveas(get(newax,'parent'),path_file,'fig')
        delete(newFigToDelete)
        return
    otherwise
        errorWindow44=errordlg(['couldn''t export because the file extension must be either: .jpg, .tif, .png, or .fig but you entered: ',path_file(end-3:end)]);
        centerFigure(hDV.hfig,errorWindow44);
        return
end

% find resolution from menu
checkedFileBinary2=strcmpi(get(hDV.uiMenuHandles.ExportFigMenu.resolutionSelectorMenus,'checked'),'on'); % find which menu option selected as on
fileResolutions=get(hDV.uiMenuHandles.ExportFigMenu.resolutionSelectorMenus,'label'); %jpeg, tiff, etc
resolutionToUse=fileResolutions(checkedFileBinary2);% find single string selected

% this line saves figure
% example: print(newFigToDelete,path_file,'-djpeg',-opengl,'-r600') % save file
%['-d',extensionToUse{1}]
print(newFigToDelete,path_file,rendererString,extensionToUse,['-r',resolutionToUse{1}]) % save file


delete(newFigToDelete)

end

% change check of selected resolution
function exportFigureTypeSettings(src,~)
% get(src,'label')
allChildrenOfParent=get(get(src,'parent'),'children'); % get all children of parent
binarySelected=allChildrenOfParent==src; % find handle of selected
set(allChildrenOfParent(binarySelected),'checked','on') % check selected
set(allChildrenOfParent(~binarySelected),'checked','off') % uncheck others
end



