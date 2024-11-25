% Stanford 2016
% Rall Walsh
%

%data entry GUI of fault data

function valsout = dataentryFaults(hDV,title,txt,vals,i,noshow)

if nargin==5
    noshow=0; %entries which we dont want to show
end

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
f = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .4*hpos(3) .7*hpos(4)];
set(f,'position',pos) ;
set(f,'MenuBar','none');
set(f,'Name',title);
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'PaperPositionMode','auto');
set(f,'closerequestfcn',@crf) ;
centerFigure(hfig,f);
set(f,'visible','on') ; 


dropDownMenu1 = uicontrol('Style', 'popup',...
    'String', {'Random Faults','Enter Faults Manually','load .CSV file'},...
    'Position', [0.3 .92 .4 .05],...
    'Callback', {@dropDownMenuFaults,hpos });

valtxtH = zeros(length(vals),1) ;

% Controls panels
% make random faults pannel
randomFaultsPannel = uipanel('parent',f,'backgroundcolor','white',    ...
    'visible','on','bordertype','none','Position',[0 0.2 1 .7]);


% make random faults pannel
pos = [0.05 0 .45 .05] ;  c2=[.6 0 -.15 0] ;
m=1 ;
for k=1:1:length(txt)
    if ~any(k==noshow)
        uicontrol('parent',randomFaultsPannel,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        
        valtxtH(k)= uicontrol('parent',randomFaultsPannel,'style','edit','string',num2str(vals(k)),...
            'pos',c2+pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        m=m+1 ;
    else
        valtxtH(k)= uicontrol('parent',randomFaultsPannel,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
end

% Controls panels
% make load Faults Pannel
loadFaultsPanel = uipanel('parent',f,'backgroundcolor','white',    ...
    'visible','off','bordertype','none','Position',[0 0.2 1 .7]);

%load faults button
loadFaultsButton = uicontrol('parent',loadFaultsPanel,'Style','push','String','Load File',...
    'pos',[.35,.3,.3,.2],'fontsize',23,'fontunits','normalized',...
    'callback',{@callbackloadfaultsfile});


uicontrol('Style','text','String',{'Load a .csv file with the following columns: X1 [km], X2 [km], Y1 [km], Y2 [km], strike [deg], dip [deg], friction coeff mu, with no header.'},'parent',...
    loadFaultsPanel,'backgroundcolor',[1 1 1],'position',[0.05,0.5,0.9,0.35]);


% type faults table
typeFaultsPannel = uipanel('parent',f,'backgroundcolor','white',    ...
    'visible','off','Position',[0 0.2 1 .7]);

if ~isfield(hDV.data.fault,'file')
    indat = cell(hDV.data.NFAULTSMAX,7);
else
    indat = cell(hDV.data.NFAULTSMAX,7);
    vals = num2cell(hDV.data.fault.file); [M,N]=size(vals) ; 
    indat(1:M,1:N)=vals; 
    
end
t = uitable('Parent', typeFaultsPannel,'units','Normalized' ,'fontunits','normalized','fontsize',.04, 'Position', [0 0 1 1], 'Data', indat);
set(t, 'ColumnName' , {'<html><h2>  &nbsp X1 [km] &nbsp   </h2></html>', '<html><h2> &nbsp X2 [km] &nbsp</h2></html>', '<html><h2> &nbsp Y1 [km]  &nbsp </h2></html>', '<html><h2>  &nbsp Y2 [km] &nbsp</h2></html>',...
    '<html><h2> &nbsp Strike [deg] &nbsp </h2></html>' , '<html><h2>  &nbsp Dip [deg]  &nbsp </h2></html>' , '<html><h2>  &nbsp Friction Coeff mu &nbsp </h2></html>'});
set(t, 'ColumnWidth', {'auto', 'auto', 'auto', 'auto', 'auto' , 'auto','auto'});
set(t, 'ColumnEditable' , [true true true true true true true]);

%OK button
uicontrol('parent',f,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH});

valsout=0;

    function crf(~,~)
        hDV.tempvec=[];
        delete(gcf);
        modal(hfig,'on');
    end


% read in csv file
    function callbackloadfaultsfile(src,~)
        
        [fileNameToLoad,pathname]=uigetfile('*.csv');
        disp(fileNameToLoad)
        readInData=csvread([pathname,fileNameToLoad]);
        set(src,'UserData',readInData) ; 
    end

    function [valtxtH]=dropDownMenuFaults(~,~,~)  % based on dropdown menu
        %      dropDownMenu
        %      source
        % callbackdata
        %         strinInputMethodSelected = get(dropDownMenu1,'String');
        index_selected = get(dropDownMenu1,'Value');
        
        switch index_selected
            case 1 % Random Faults
                %                 pannel1 on
                % pannel 2,3 off visibility
                set(loadFaultsPanel,'visible','off')
                set(typeFaultsPannel,'visible','off')
                set(randomFaultsPannel,'visible','on')
                hpos = get(hfig,'position') ;
                hpos2 = get(f,'position') ;
                pos = [(hpos2(1)) (hpos2(2)+hpos2(4)-.7*hpos(4)) .4*hpos(3) .7*hpos(4)];             
                set(f,'position',pos) ;
            case 2 % Enter Faults Manually
                disp('enter faults')
                set(loadFaultsPanel,'visible','off')
                set(randomFaultsPannel,'visible','off')
                set(typeFaultsPannel,'visible','on')
                hpos = get(hfig,'position') ; % original position
                hpos2 = get(f,'position') ; % current position of this figure
                pos = [(hpos2(1)) (hpos2(2)+hpos2(4)-.7*hpos(4)) .4*hpos(3) .7*hpos(4)];   
                set(f,'position',pos) ;
            case 3 % load .CSV file
                set(randomFaultsPannel,'visible','off')
                set(typeFaultsPannel,'visible','off')
                set(loadFaultsPanel,'visible','on')
                hpos = get(hfig,'position') ;
                %hpos2 = get(f,'position') ; 
                %pos = [hpos2(1) (hpos2(2)+hpos2(4).*.5) hpos2(3) hpos2(4).*.5];
                %set(f,'position',pos) ;
                disp('load')
        end
        
        valtxtH=[];
    end



%this must tie back to callbackdatabuts.m
    function butcall(~,~,hDV,htxt)
        inputMode = get(dropDownMenu1,'Value'); 
        vec=zeros(length(htxt),1);
        for n=1:1:length(htxt)
            vec(n)=str2double(get(htxt(n),'string'));
        end
        switch inputMode
            case 1
                % random faults
                hDV.data.fault.vals=vec;
            case 2
                hDV.data.fault.file=cell2mat(get(t,'data'));
            case 3
                % faults for file
                hDV.data.fault.file=get(loadFaultsButton,'UserData') ; 
        end
        
        hDV.data.fault.intype = inputMode;
        
        if inputMode==2 || inputMode==3
            [~,cols]=size(hDV.data.fault.file) ;
            if any(isnan(hDV.data.fault.file(:))) || cols~=7
                ed=errordlg('Check Data Entry for Faults','Faults Entry','modal') ;
                waitfor(ed);
            else
                delete(gcf);
                modal(hfig,'on');
            end
        end
        
        if inputMode==1
            delete(gcf);
            modal(hfig,'on');
        end
        
    end


end




