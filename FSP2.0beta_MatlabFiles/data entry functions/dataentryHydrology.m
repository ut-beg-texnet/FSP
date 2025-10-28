% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Standard data entry GUI

function valsout = dataentryHydrology(hDV,title,txt,vals,i,noshow)

if nargin==5
    noshow=0; %entries which we dont want to show
end

if hDV.plotdata.printFunctionName
    disp(['running dataentryHydrology'])
end
hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
loadHydrologyFigure = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .5*hpos(3) .8*hpos(4)];
set(loadHydrologyFigure,'position',pos) ;
set(loadHydrologyFigure,'MenuBar','none');
set(loadHydrologyFigure,'Name',title);
set(loadHydrologyFigure,'NumberTitle','off');
set(loadHydrologyFigure,'DefaultUicontrolUnits','normalized');
set(loadHydrologyFigure,'DefaultUicontrolFontsize',14);
set(loadHydrologyFigure,'PaperPositionMode','auto');
set(loadHydrologyFigure,'closerequestfcn',@crf) ;
centerFigure(hfig,loadHydrologyFigure);
set(loadHydrologyFigure,'visible','on') ;
pos = [0.05 0 .45 .05] ;  c2=[.6 0 -.15 0] ;
valtxtH = zeros(length(vals),1) ;



radioButtonGroupHydrology=uibuttongroup('Position', [0 .85 .8 .18],'backgroundcolor',[1 1 1],...
    'SelectionChangeFcn', {@dropDownCallbackHydrology},'parent',loadHydrologyFigure,'fontsize',hDV.ftsz,...
    'BorderType','none','fontunits','normalized');

radioButton1=uicontrol(radioButtonGroupHydrology,'Style','radiobutton',...
    'Position', [0.01 .45 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Enter Hydrologic Parameters','fontunits','normalized',...
    'ToolTipString','Select to randomly generate faults');

radioButton2=uicontrol(radioButtonGroupHydrology,'Style','radiobutton',...
    'Position', [0.01 .0 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Load External Hydrologic Model','fontunits','normalized',...
    'ToolTipString','Select to load a .csv file');


% % dropdown menu, decide to load file, or enter wells manually
% dropDownMenuHydrology = uicontrol('Style', 'popup',...
%     'String', {'Enter Hydrologic Parameters' ,  'Load External Hydrology Model'},...
%     'Position', [0.1 .8 .8 .15],'fontsize',18,'fontunits','normalized',...
%     'Callback', @dropDownCallbackHydrology);

if  isfield(hDV.data.reservoir,'importHydrology') % if has been selected already
    if hDV.data.reservoir.importHydrology==1
        set(radioButtonGroupHydrology,'SelectedObject',radioButton2);
    end
end

enterHydroParametersPanel = uipanel('parent',loadHydrologyFigure,'backgroundcolor','white',    ...
    'visible','off','bordertype','none','Position',[0 0.2 1 .7]);

m=1 ;
for k=1:1:length(txt)
    if ~any(k==noshow)
        uicontrol('parent',enterHydroParametersPanel,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        
        valtxtH(k)= uicontrol('parent',enterHydroParametersPanel,'style','edit','string',num2str(vals(k)),...
            'pos',c2+pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        m=m+1 ;
    else
        valtxtH(k)= uicontrol('parent',enterHydroParametersPanel,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
end

%%%%%%%%%%%%%%%%%%
% load hydrology xyzpt
%%%%%%%%%%%%%%%%
% panel
loadExternalHydroPanel = uipanel('parent',loadHydrologyFigure,'backgroundcolor','white',    ...
    'visible','off','bordertype','none','Position',[0 0.2 1 .7]);

% load wells button
uicontrol('parent',loadExternalHydroPanel,'Style','push','String','Load .csv File',...
    'pos',[0.6,0.85,.2,.1],'fontsize',14,'fontunits','normalized',...
    'callback',{@callbackloadadvancedHuydrologyfile});

%text
uicontrol('Style','text','String',{'Number of header lines:'} ... % ;' API#, X [km], Y [km], Year, month (1-12), Injection volume [bbl/month], Pressure [psi], Well name, and well type' }...
    ,'parent',loadExternalHydroPanel,'fontsize',16,'fontunits','normalized','HorizontalAlignment','right'... % ,'HorizontalAlignment','left'
    ,'backgroundcolor','white','position',[0.02,0.86,0.45,0.06]);

% number of headers in file to load
uicontrol('parent',loadExternalHydroPanel,'style','edit','string',num2str(hDV.data.reservoir.loadedHeaderLines),...
    'pos',[.51 .87 .05 .06],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'callback',{@nHydrologyHeaderscall});

% make uitable panel
uitableHydrologyData = uitable('parent',loadExternalHydroPanel,'units','normalized','Position',[0.01 0.1 0.98 .72],'ColumnName', {'East km';'North km';['Change in PSI at ',num2str(round(hDV.data.stress.vals(5))),'ft'];'Year (Jan 1, 12:01 AM)'},...
    'ColumnEditable', [false(4,1)]',...  'data',hDV.data.realWellData.stringsWellDataAdvanced,
    'tag','uitableWellDataTag1','CellEditCallback',{@editedUITableWellData});

if  isfield(hDV.data.reservoir,'stringsImportedHydrology') % if has been selected already, set dropdown menu value
    set(uitableHydrologyData,'data', hDV.data.reservoir.stringsImportedHydrology)
end

%help format well data button
uicontrol('parent',loadExternalHydroPanel,'style','push','string','File Format Help',...
    'position',[0.0 0.0  0.4 .1],'fontsize',14,'fontunits','normalized','callback',{@callbuthelpformatHydrology,loadHydrologyFigure});

%resize column name header font size
hs = '<html><h2>'; %html start
he = '</h2></html>'; %html end
cn = get(uitableHydrologyData,'ColumnName'); %original
cnh = cellfun(@(x)[hs x he],cn,'uni',false); %with html
set(uitableHydrologyData,'ColumnName',cnh) %apply


%%%%%%%%%%%




dropDownCallbackHydrology([],[])

%OK button
uicontrol('parent',loadHydrologyFigure,'style','pushbutton','string','OK',...
    'pos',[.25 .01 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,i,loadHydrologyFigure});

valsout=0;

    function crf(~,~)
        hDV.tempvec=[];
        delete(gcf);
        modal(hfig,'on');
    end

%this must tie back to callbackdatabuts.m
    function butcall(~,~,hDV,htxt,i,loadHydrologyFigure)
        vec=zeros(length(htxt),1);
        for n=1:1:length(htxt)
            vec(n)=str2double(get(htxt(n),'string'));
        end
        
        %         hDV.data.reservoir.importHydrology=get(dropDownMenuHydrology,'val')-1;
        if  get(radioButtonGroupHydrology,'SelectedObject')==radioButton1;
            hDV.data.reservoir.importHydrology=0;
        else
            hDV.data.reservoir.importHydrology=1;
        end
        
        
        
        resetButtons=0;
        switch   hDV.data.reservoir.importHydrology
            case 0 % doing hydrology calculation with wells
                if ~all(hDV.data.reservoir.vals==vec) ; resetButtons=1;end
                hDV.data.reservoir.vals=vec;
            case 1 % imported MODFLOW model
                if isempty(get(uitableHydrologyData,'data'))
                    msgWindow1=msgbox({cat(2,'You didn''t enter a hydrologic model. Switch back to select entering hydrologic parameters if you don''t have an external model to bring in. ')}, 'Loading Hydrologic Model','warn');
                    centerFigure(hDV.hfig,msgWindow1);
                    hDV.data.reservoir.importHydrology=0; %
                else % get data
                    SpreadsheetStrings2HydrologyData(get(uitableHydrologyData,'data'),hDV); % function that reads uitable data and saves it to structure
                    resetButtons=1;
                     hDV.data.probHydrology.probabilistic1VsDeterministic2=2;
                end
        end
        
        if resetButtons==1
            % then some data was changed
            % Set Calculate and geomechanics button red
            % don't need to reset prob geomechanics button
            resetButtonsRed(hDV,[1,3,4])
        end
        
        %                         if ~isempty(hDV.data.realWellData.stringsWellDataAdvanced) % don't call function if no data is loaded in
        %                     SpreadsheetStrings2WellData(get(uitableWellData,'data'),hDV); % function that reads uitable data and saves it to structure
        %                         end
        
        
        delete(loadHydrologyFigure);
        modal(hfig,'on');
    end



    function dropDownCallbackHydrology(src,~)  % based on dropdown menu, toggle visibility

        %         index_selected = get(dropDownMenuHydrology,'Value');
        %         switch index_selected
        %             case 1 % Enter hydrologic parameters
        if get(radioButtonGroupHydrology,'selectedobject') ==radioButton1
            %                 disp('enter wells manually chosen')
            set(enterHydroParametersPanel,'visible','on')
            set(loadExternalHydroPanel,'visible','off')
            hDV.data.reservoir.importHydrology=0;
            
            %             case 2 % csv import x,y,p,t
        else
            %                 disp('enter wells csv chosen')
            set(enterHydroParametersPanel,'visible','off')
            set(loadExternalHydroPanel,'visible','on')
            hDV.data.reservoir.importHydrology=1;
        end
    end

% read in csv file with monthly well data
    function callbackloadadvancedHuydrologyfile(~,~) % read well data csv
        [fileNameToLoad,pathname]=uigetfile('*.csv');
        if any(fileNameToLoad~=0)
            fid44=fopen([pathname,fileNameToLoad]); % open file
            readInData=textscan(fid44,'%s%s%s%s','Delimiter',',','Headerlines',hDV.data.reservoir.loadedHeaderLines,'multipledelimsasone',true); % read all data in as text strings
            fclose(fid44) ; % close file
            %
            hDV.data.reservoir.stringsImportedHydrology=[readInData{:,:}];
            
            [hDV,flagHydrology]=checkHydrologyData(hDV); % check data, throw flags
            if flagHydrology~=1 % if no issues with data, then proceed to load it
                set(uitableHydrologyData,'data', hDV.data.reservoir.stringsImportedHydrology)
                set(uitableHydrologyData,'ColumnEditable',[false(size( hDV.data.reservoir.stringsImportedHydrology,2),1)]')
            end
        end
    end

% callback reset headers
    function nHydrologyHeaderscall(src,~)
        val = str2double(get(src,'string')) ;
        val=floor(val);
        if val<0
            val=0;
        end
        set(src,'string',num2str(val));
        hDV.data.reservoir.loadedHeaderLines=val;
    end
% help formatting well data image function
    function callbuthelpformatHydrology(~,~,hfig)
        modal(hfig,'off');
        hpos = get(hfig,'position') ;
        helpImageFigure = figure('color',[1 1 1],'units','pixels','visible','on');
        pos = [0 0 .6*hpos(3) .8*hpos(4)];
        imshow('HydrologyDataLoadHelp.png','border','tight');
        set(helpImageFigure,'position',pos) ;
        set(helpImageFigure,'MenuBar','none');
        set(helpImageFigure,'NumberTitle','off');
        set(helpImageFigure,'DefaultUicontrolUnits','normalized');
        set(helpImageFigure,'DefaultUicontrolFontsize',14);
        set(helpImageFigure,'closerequestfcn',@crf) ;
        centerFigure(hfig,helpImageFigure);
        %set(f,'visible','on') ;
        function crf(~,~)
            delete(helpImageFigure);
            modal(hfig,'on');
        end
        
    end



end




