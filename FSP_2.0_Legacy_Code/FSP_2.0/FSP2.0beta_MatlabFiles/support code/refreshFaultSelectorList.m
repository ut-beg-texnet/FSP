% Update listbox of faults selected
% take in current tab, HDV
% get: number of faults, color based on tab and appropriate colormap
% update maximum number that can be selected
% update hDV.curfault value

% but do not invoke callback of listbox function that is being updated, because that calls
% refreshplotdata, which this function is called from.

% By Rall Walsh
% Stanford

function refreshFaultSelectorList(hDV,currentTabName)
% setup for coding, and get current tab if necessary
if nargin == 0
    hDV=evalin('base', 'hSV'); % get hDV to test?
    currentTabName=hDV.currtab.name;
elseif nargin ==1
    currentTabName=hDV.currtab.name;
end

% get numbers from hDV
menuHandle=hDV.plotdata.ListboxFaultSelector; % handle of uicontrol object
nfaults=hDV.data.fault.vals(1); % number of faults
curfault=hDV.plotdata.curfault ;% current fault(s)


% if ppfail or FSP not yet calculated, then don't color
%
if any(strcmp({'GEOMECHANICS','PROB. GEOMECH','HYDROLOGY'},currentTabName)) && ~isfield(hDV.plotdata.results.outs,'ppfail')
    currentTabName='MODEL INPUTS'; % this will only prevent coloring below; will not affect tabname elsewhere
elseif any(strcmp({'INTEGRATED','PROB. HYDRO'},currentTabName)) && ~isfield(hDV.plotdata.pint,'fsp')
    currentTabName='MODEL INPUTS'; % this will only prevent coloring below; will not affect tabname elsewhere
end

%
% Color faults in fault list on the left of main window
%

% Initially don't color faults in the list; just setup fault number text strings
%
xstr = cell(nfaults+1,1) ; xstr{1}='All Faults';
for k=1:1:nfaults
    xstr{k+1}=['Fault #',num2str(k)];
end
set(menuHandle,'string',xstr)

% Now color faults only if calculate button is NOT red
%
if ~all(get(hDV.bCalc,'backgroundcolor')==hDV.colors.red)
    switch currentTabName   
        case {'GEOMECHANICS','PROB. GEOMECH','HYDROLOGY'}
            faultFSPToColorby=hDV.plotdata.results.outs.ppfail; % Use ppfail for color instead of FSP here

            minint=hDV.plotdata.mincfc; % color bar limits
            maxint=hDV.plotdata.maxcfc;
            probHydDropdownTxtClr(menuHandle,hDV,faultFSPToColorby,minint,maxint,flipud(hDV.cmapGYR),0)

        case {'INTEGRATED','PROB. HYDRO'} 
            faultFSPToColorby=hDV.plotdata.pint.fsp; % FSP value

            minint=hDV.plotdata.minint; % color bar limits
            maxint=hDV.plotdata.maxint;
            probHydDropdownTxtClr(menuHandle,hDV,faultFSPToColorby,minint,maxint,hDV.cmapGYR,1) % set colors of menu

        otherwise  % case 'MODEL INPUTS'

    end
end
% update menu values
set(menuHandle,'max',nfaults+1,'value',curfault+1);

end

