

% By Rall
% Make input data uimenu
% Stanford 2017

function makeInputDataUiMenu(hDV)


hDV.uiMenuHandles.DataInputMenu.Tophandle=uimenu(hDV.hfig,'Label','Data Inputs');

% top tabs
hDV.uiMenuHandles.DataInputMenu.dbutnames = {'Stress Data','Hydrology Data','Injection Wells','Fault Data','Advanced'};
%             dataInputButtonTooltipStrings={'Enter Stress Field data','Enter Hydrologic data - Porosity, permeability, etc.',...
%                 'Enter Injection well data, either constant rates, or monthly injection rates','Randomly generate faults, or enter fault strike, dip, location, etc.'...
%                 'additional parameters that can be edited'};
%             pos=[.05 .9 .9 .05] ; del=[0 .06 0 0] ;


for i = 1:1:length(hDV.uiMenuHandles.DataInputMenu.dbutnames)
    hDV.uiMenuHandles.DataInputMenu.dButs(i) = uimenu('parent',hDV.uiMenuHandles.DataInputMenu.Tophandle,...
        'Label',hDV.uiMenuHandles.DataInputMenu.dbutnames{i},...
        'callback',{@callbackdatabuts,i,hDV});
end
end