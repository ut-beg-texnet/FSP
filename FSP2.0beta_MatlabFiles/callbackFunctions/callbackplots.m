function callbackplots(~,~,hDV,num,name)

% Calculation tab: either prescribed on input or remains the same
if ( ~isempty(num) &&  ~isempty(name) )
    hDV.currtab.number = num;
    hDV.currtab.name   = name;
else
    num  = hDV.currtab.number;
    name = hDV.currtab.name;
end

if hDV.plotdata.printFunctionName
    disp(['running callbackplots ',name])
end

% Heights of buttons at top of screen
HB    = hDV.bMain.bTopHeights.HB;
HBdel = hDV.bMain.bTopHeights.HBdel;

% Set tab sizes, colors and fonts
for i = 1:1:length(hDV.bMain.bTop)
    pos = get(hDV.bMain.bTop(i),'position');
    if ( i ~= num )
        pos(4) = HB;
        col    = hDV.colors.darkgrey;
        set(hDV.bMain.bTop(i),'position',pos,'foregroundcolor',col,'fontweight','bold');
        set(hDV.pMain(i),'visible','off');
    else
        pos(4) = HB + HBdel;
        col    = hDV.colors.black;
        set(hDV.bMain.bTop(i),'position',pos,'foregroundcolor',col,'fontweight','demi');
        set(hDV.pMain(i),'visible','on'); %turn on the selected panel
        calcengine(hDV,name); %calculation engine runs for each panel
        refreshplotdata(hDV,name); %plot data for selected panel refreshed
    end
    if i<=num+1 % only enable button after previous has been clicked
        set(hDV.bMain.bTop(i),'Enable','on');
    end
end

%enable top buttons
% set(hDV.bMain.bTop(:),'Enable','on');
set(hDV.pMain(end),'visible','off'); %this is the blank panel that is turned off; defined in surfaceviz class

end