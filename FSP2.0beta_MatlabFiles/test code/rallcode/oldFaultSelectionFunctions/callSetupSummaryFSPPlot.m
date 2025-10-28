% callback to make summary plot
% by Rall Walsh
% added 0.98.6
% to do: export data button
% allow selection of deterministic hydrology
% 


% replaced by callRunSummaryPlot in 0.98.9

  function callSetupSummaryFSPPlot(~,~,hDV)
title='Generate Summary Plot (can take a while with many years)';
txt={'Start Year';'End Year';'Slip Potential Y limit (0 to 1)';'0 for deterministic, 1 for probabilistic hydrology'};
vals=[get(hDV.hdsldr(1),'min');get(hDV.hdsldr(1),'max');1;0];
stringsTooltip={'Whole year to start making summary plot';'Whole year to end making summary plot';...
    'Vertical axis limit of Fault FSP through time plot';'Choose deterministic or probabilistic hydrology'};

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
figureEntryHandle159 = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .4*hpos(3) .5*hpos(4)];
set(figureEntryHandle159,'position',pos) ;
set(figureEntryHandle159,'MenuBar','none');
set(figureEntryHandle159,'Name',title);
set(figureEntryHandle159,'NumberTitle','off');
set(figureEntryHandle159,'DefaultUicontrolUnits','normalized');
set(figureEntryHandle159,'DefaultUicontrolFontsize',14);
set(figureEntryHandle159,'PaperPositionMode','auto');
set(figureEntryHandle159,'closerequestfcn',@crf) ;
centerFigure(hfig,figureEntryHandle159);
set(figureEntryHandle159,'visible','on') ;
pos = [0.0 0 .6 .18] ;  c2=[.65 0.1 -.3 -.1] ;
% valtxtH = zeros(length(vals),1) ;

m=1 ;

selectionsPanel= uipanel('parent',figureEntryHandle159,'backgroundcolor','white',    ...
    'visible','on','bordertype','none','Position',[0 0.15 .7 .85]);
for k=1:1:3 %length(txt) % this allows third row - deterministic or probabilistic designation
    
    uicontrol('parent',selectionsPanel,'style','text','string',txt{k},...
        'pos',pos+[0 .95-m*.2 0 0],'fontsize',12,'fontunits','normalized',...
        'tooltipstring',stringsTooltip{k},'backgroundcolor',[1 1 1]);
    
    valtxtH(k)= uicontrol('parent',selectionsPanel,'style','edit','string',num2str(vals(k)),...
        'pos',c2+pos+[0 .95-m*.2 0 0],'fontsize',12,'fontunits','normalized',...
        'tooltipstring',stringsTooltip{k},'backgroundcolor',[1 1 1]);
    m=m+1 ;
end

nfaults=hDV.data.fault.vals(1);
% make multiple fault selections
%setup fault number text
xstr = cell(nfaults+1,1) ; xstr{1}='All' ;
for k=1:1:nfaults
    xstr{k+1}=['Fault #',num2str(k)];
end

uicontrol('parent',figureEntryHandle159,'style','text','string','Select Fault(s)','position',[.7 .9 .3 .1],...
    'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1]);

uiListboxHandle=uicontrol('style','listbox','parent',figureEntryHandle159,'string',xstr,'position',[.7 .0 .3  .9],...
    'fontsize',12,'fontunits','normalized','value',1,'backgroundcolor',[1 1 1],...
    'min',1,'max',nfaults,...
    'TooltipString','hold control to select multiple faults, hold shift to select multiple faults in a row');



%OK button
uicontrol('parent',figureEntryHandle159,'style','pushbutton','string','Go!',...
    'pos',[.25 .05 .45 .1],'fontsize',12,'fontunits','normalized',...
    'TooltipString','Press to cycle over calculation one year at a time and make summary plot',...
    'callback',{@butcall,hDV,valtxtH,figureEntryHandle159,uiListboxHandle});

    function crf(~,~)
        delete(gcf);
        modal(hfig,'on');
    end

end

function butcall(~,~,hDV,valtxtH,figureEntryHandle159,uiListboxHandle)

valsout=zeros(length(valtxtH),1); % get field values
for n=1:1:length(valtxtH)
    valsout(n)=str2double(get(valtxtH(n),'string'));
end
startyear=get(hDV.hdsldr(1),'min');
endyear=get(hDV.hdsldr(1),'max'); % this is hardwired in slider bar

valuesSelected=get(uiListboxHandle,'value')-1;
nfaults=hDV.data.fault.vals(1);

if any (valuesSelected==0) % if all was selected (at all)
    faultsSelected=1:nfaults; 
else % just certain faults
    faultsSelected=valuesSelected;
end

if valsout(1)>valsout(2) || valsout(1)< startyear || valsout(2) > endyear
    errorWindow1=errordlg(cat(2,' check the years you entered :(  '));
    centerFigure(hDV.hfig,errorWindow1);
    return % don't save
end

if valsout(3)>1 || valsout(3)<=0
    errorWindow1=errordlg(cat(2,' check the FSP axis limit value you entered, it should be between 0 and 1, but it''s:  ',num2str(valsout(3))));
    centerFigure(hDV.hfig,errorWindow1);
        return % don't save
end

YearsToRun=valsout(1):1:valsout(2);
SavedFSPmatrix=zeros(length(hDV.plotdata.pint.fsp),length(YearsToRun));

% close data entry handle
delete(figureEntryHandle159);
modal(hDV.hfig,'on');


handleWaitbar1=waitbar(0.01,'Calculating FSP January 1 each year','CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
centerFigure(hDV.hfig,handleWaitbar1)
setappdata(handleWaitbar1,'canceling',0)
runCancelled=0;
for kjj449=1:length(YearsToRun)
    
    waitbar(kjj449./length(YearsToRun),handleWaitbar1,['Calculating FSP in ',num2str(YearsToRun(kjj449))]);
    pause(1) % 1 second
    set(hDV.hdsldr(3),'val',YearsToRun(kjj449)) % update waitbar value
    sldrcall(hDV.hdsldr(3),[],hDV) % run calculation
    if getappdata(handleWaitbar1,'canceling') % if cancel button is hit, stop
        disp('Summary plot generation cancelled')
        runCancelled=1;
        break
    end
    SavedFSPmatrix(:,kjj449)=hDV.plotdata.pint.fsp;  % save value
end
delete(handleWaitbar1)  % DELETE the waitbar; don't try to CLOSE it.

% call make plot function
if ~runCancelled
    PlotSummaryPlot(SavedFSPmatrix,YearsToRun,hDV,faultsSelected,valsout(3))
end

end


function PlotSummaryPlot(SavedFSPmatrix,YearsToRun,hDV,faultsSelected,yLimit22)
% take data that was calculated and make plot
hfig = hDV.hfig;

% hpos = get(hfig,'position') ;
% figureEntryHandle157 = figure('tag','Summary FSP through time','color',[1 1 1],'units','pixels','visible','off');
% pos = [0 0 .4*hpos(3) .5*hpos(4)];
% set(figureEntryHandle157,'position',pos) ;
% set(figureEntryHandle157,'MenuBar','none');
% set(figureEntryHandle157,'Name','Fault Slip Potential Through Time');
% set(figureEntryHandle157,'NumberTitle','off');
% set(figureEntryHandle157,'DefaultUicontrolUnits','normalized');
% set(figureEntryHandle157,'DefaultUicontrolFontsize',14);
% set(figureEntryHandle157,'PaperPositionMode','auto');
% set(figureEntryHandle157,'closerequestfcn',{@crf22,figureEntryHandle157}) ;
% centerFigure(hfig,figureEntryHandle157);
% set(figureEntryHandle157,'visible','on') ;

visibility1='on';

%plot fsp with time on each fault
% pint.ax4=axes('parent',figureEntryHandle157,'position',[.1  .2  .8   .72],'visible',visibility1) ; % ,'fontsize',hDV.ftsz
pint=hDV.plotdata.pint;

title(pint.ax4,'Fault Slip Potential Through Time');
hold(pint.ax4,'on') ;
xlabel(pint.ax4,'Time [years]','fontsize',hDV.ftsz) ; ylabel(pint.ax4,'Fault Slip Potential','fontsize',hDV.ftsz) ;
set(pint.ax4,'ylim',[0,yLimit22]);
startAxis=min(YearsToRun)-0.5;
endAxis=max(YearsToRun)+0.5;
set(pint.ax4,'xlim',[startAxis,endAxis]);

% background color patches
set(hDV.plotdata.pffot.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
colorBarPatches=zeros(length(hDV.cmapInfo.patchesInColorbar),1);
patchYs=linspace(0,1, hDV.cmapInfo.patchesInColorbar+1); % Y coordinates of bounds of patches
PatchXCoord= [startAxis,endAxis,endAxis,startAxis];% X coordinates of patches
for j=1:hDV.cmapInfo.patchesInColorbar % cycle over each patch
    PatchYCoord = [patchYs(j),patchYs(j),patchYs(j+1),patchYs(j+1)] ;
    colorBarPatches(j,1)=patch(PatchXCoord,PatchYCoord,[0 0 0],'parent',pint.ax4,'edgecolor','none','visible','off') ;
    
end

%get the numbers from GUI
minint = str2double(get(hDV.plotdata.pint.cmintxt,'string'));
maxint = str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));
%check and reset if needed
if isnan(minint) || isnan(maxint) || minint>=maxint || minint<0 || maxint>1 % if variables don't exist correctly, do nothing for now
    
else
    % make color patches
    cv=linspace(0,1,hDV.cmapInfo.patchesInColorbar);
    %edit colorbar range
    for j=1:hDV.cmapInfo.patchesInColorbar
        
        cl = getcolor(hDV.cmapGYR,cv(j),minint,maxint) ;
        set(colorBarPatches(j), 'facecolor', cl,'visible','on');
        
    end
end


ticks=get(pint.ax4,'xtick');
ticks=ticks(ticks==round(ticks)); % remove non whole years
set(pint.ax4,'xtick',ticks);


pint.FaultFSPCurvesThruTime=zeros( hDV.data.fault.vals(1),1);

% make curves of FSP vs time for each fault
for j=faultsSelected
    labelYcoord=[];
    labelXcoord=[];
    % this is setup code
    pint.FaultFSPCurvesThruTime(j) = plot(0,0,'parent',pint.ax4,...
        'linewidth',2,'HandleVisibility','off','color','b','visible',visibility1,...
        'marker','x');
    % this is refresh code:
    set(pint.FaultFSPCurvesThruTime(j),'xdata',YearsToRun,'ydata',SavedFSPmatrix(j,:))
    
    % decide which points to label
    for kindx445=1:size(SavedFSPmatrix,2) % cycle over year index
        [uniqueValues,faultsWithUniqueValues,~]=unique(SavedFSPmatrix(:,kindx445),'first'); % find first instance of each FSP value
        if any(faultsWithUniqueValues==j) && uniqueValues(faultsWithUniqueValues==j)>0.0011 % then we want a label at this location and year
            % also make a requirement that value is greater than 1%?
            labelYcoord(end+1,1)=uniqueValues(faultsWithUniqueValues==j); % save label coordinates
            labelXcoord(end+1,1)=YearsToRun(kindx445);
        end
    end
    
    
    % downsample if lots of years, try to stagger too
    if length(labelXcoord)>=41 && ~isempty(labelXcoord)
        labelXcoord=labelXcoord(mod(j,5)+1:8:end,1);
        labelYcoord=labelYcoord(mod(j,5)+1:8:end,1);
    elseif length(labelXcoord)>=21 && ~isempty(labelXcoord)
        labelXcoord=labelXcoord(mod(j,3)+1:4:end,1);
        labelYcoord=labelYcoord(mod(j,3)+1:4:end,1);
    elseif length(labelXcoord)>=8 && ~isempty(labelXcoord)
        labelXcoord=labelXcoord(mod(j,2)+1:2:end,1);
        labelYcoord=labelYcoord(mod(j,2)+1:2:end,1);
    end
    
    text(labelXcoord,labelYcoord,num2str(j),'verticalalignment','bottom','color','b')
    
%     
end
%time bars on the right plots
grid(pint.ax4,'on');
% %button to export CSV of FSP vs time curves% needs to be adapted for this
% plot
 uicontrol('parent',figureEntryHandle157,'style','push','string',{'Export data'},'position',[.0   .0   .2  .08],...
    'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callExportFSPThruoghTime_CSV,hDV,pint,faultsSelected},'visible',visibility1);

set(hDV.plotdata.pprob.summaryPlotButton,'backgroundcolor',hDV.colors.green)

    function crf22(~,~,figureEntryHandle157)
        delete(figureEntryHandle157);
    end
end

function callExportFSPThruoghTime_CSV(~,~,hDV,pint,faultsSelected)

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
    
    fprintf(fid22, 'fault data:,,1 fault per line,, and 1 year per column ,,,Fault Slip Potential (units of probability, in same column as corresponding date). Date is January 1st at 12:01 AM in corresponding year\n');
    fclose(fid22);% close
    
    % get year data fom the first fault that's calculated
    k=pint.FaultFSPCurvesThruTime(pint.FaultFSPCurvesThruTime~=0);
    k=k(1);
    TimeData=get(k,'xdata');
    
    
    dlmwrite(path_file,TimeData,'-append','coffset',7)% write data
    fid223 = fopen(path_file, 'a'); % append
    fprintf(fid223, ['Fault Number, CenterX, CenterY, Strike(deg), Dip(Deg), length(km),Coefficient of Friction (mu),Fault Slip Potential in that year,,\n']); % headers
    fclose(fid223);
    

    for k443Idx= 1:length(faultsSelected) % cycle over faults to print data
        j443Idx=faultsSelected(k443Idx); % changed this from indexing variable when only outputting selected faults
        thisCurve=pint.FaultFSPCurvesThruTime(j443Idx) ;        %pick fault curve
        PSIData=get(thisCurve,'ydata');% this fault PSI data
        %                RGBDataThisFault=get(thisCurve,'color'); % when was saving color too
        dlmwrite([path2,fileNameToPrint],[j443Idx,hDV.data.fault.xf(j443Idx),hDV.data.fault.yf(j443Idx),hDV.data.fault.thf(j443Idx),hDV.data.fault.dipf(j443Idx),hDV.data.fault.lenf(j443Idx)...
            ,hDV.data.fault.muf(j443Idx) ,PSIData],'-append');  % ,'coffset',2),RGBDataThisFault(1),RGBDataThisFault(2),RGBDataThisFault(3) ,red (RGB 0to1),green (RGB 0to1),blue (RGB 0to1),
    end% cycle over faults to print data
end % if file
end % function



function sldrcall(src,~,hDV)
val = get(src,'value') ; val = floor(val) ;
mival = get(src,'min') ; mxval = get(src,'max') ;
if val<mival
    val=mival;
elseif val>mxval
    val=mxval;
end
set(src,'value',val) ;
%         set(yrtxt,'string',num2str(val)) ;

%update all other sliders and text
set(hDV.hdsldr(:),'value',val) ;
set(hDV.hdsldr_txt(:),'string',num2str(val)) ;

calcengine(hDV,hDV.currtab.name); %rerun calculations
refreshplotdata(hDV,hDV.currtab.name); %plot data for selected panel refreshed
end
