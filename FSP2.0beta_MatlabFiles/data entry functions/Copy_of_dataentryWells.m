% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

% Grid data entry GUI
% modified by Rall Walsh to also load csv files


function valsout = dataentryWells(hDV,title,txt,vals,~)

hfig = hDV.hfig;
modal(hfig,'off');

[M,N] = size(vals) ; %dimensions of data entry matrix; rows are individual wells
if M>8, M=8; end     %regardless of the maximum number of wells, set the max for entering wells to 8.
valtxtH = zeros(M,N) ;

hpos = get(hfig,'position') ;
loadWellsFigure = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .7*hpos(3) .9*hpos(4)];
set(loadWellsFigure,'position',pos) ;
set(loadWellsFigure,'MenuBar','none');
set(loadWellsFigure,'Name',title);
set(loadWellsFigure,'NumberTitle','off');
set(loadWellsFigure,'DefaultUicontrolUnits','normalized');
set(loadWellsFigure,'DefaultUicontrolFontsize',14);
set(loadWellsFigure,'PaperPositionMode','auto');
set(loadWellsFigure,'closerequestfcn',@crf) ;
centerFigure(hfig,loadWellsFigure);
set(loadWellsFigure,'visible','on') ;

% dropdown menu, decide to load file, or enter wells manually
dropDownMenu2 = uicontrol('Style', 'popup',...
    'String', {'Enter Wells Manually' ,  'Load Wells Complete .csv'},...
    'Position', [0.2 .8 .6 .15],'fontsize',18,'fontunits','normalized',...
    'Callback', @dropDownMenuWells);

if isfield(hDV.data.realWellData,'use')
    set(dropDownMenu2,'Value',hDV.data.realWellData.use+1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter wells one at a time panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enterWellsPanel = uipanel('parent',loadWellsFigure,'backgroundcolor','white',    ...
    'visible','off','bordertype','none','Position',[0 0.2 1 .7]);

% label wells number field
uicontrol('parent',enterWellsPanel,'style','text','string','Number of wells:',...
    'pos',[.15 .9 .4 .05],'fontsize',14,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);

% shorthand well label
pos0=[.05 0 .1 .05];
for y=1:1:M
    uicontrol('parent',enterWellsPanel,'style','text','string',num2str(y),... % char(65+y-1),...
        'pos',pos0+[0 .8-y*.08 0 0],'fontsize',12,'fontunits','normalized',...
        'backgroundcolor',[1 1 1]);
end

wd=.8/length(txt) ; pos = [0 .8 wd .05] ;
for x=1:1:N
    uicontrol('parent',enterWellsPanel,'style','text','string',txt{x},...
        'pos',pos+[.15+(x-1)*wd 0 0 0],'fontsize',12,'fontunits','normalized',...
        'backgroundcolor',[1 1 1]);
end


for m=1:1:M
    for k=1:1:N
        posoff = [0 -m*.08  0 0] ;
        valtxtH(m,k)= uicontrol('parent',enterWellsPanel,'style','edit','string',num2str( hDV.data.inject.vals(m,k)),...
            'pos',pos+[.15+(k-1)*wd 0 0 0]+posoff,'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
    end
end


% number of wells field
hNwells = uicontrol('parent',enterWellsPanel,'style','edit','string',num2str(hDV.data.nwells),...
    'pos',[.45 .9 .2 .05],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'callback',{@nwellscall,valtxtH,M});

%setup to start
nwellscall(hNwells,[],valtxtH,M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make load advanced wells Panel2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadWellsPanel2 = uipanel('parent',loadWellsFigure,'backgroundcolor','white',    ...
    'visible','off','bordertype','none','Position',[0 0.15 1 .75],'visible','off');

% load wells button
uicontrol('parent',loadWellsPanel2,'Style','push','String','Load .csv File',...
    'pos',[0.6,0.85,.2,.1],'fontsize',14,'fontunits','normalized',...
    'callback',{@callbackloadadvancedwellsfile});

%text
uicontrol('Style','text','String',{'Number of file header lines:'} ... % ;' API#, X [km], Y [km], Year, month (1-12), Injection volume [bbl/month], Pressure [psi], Well name, and well type' }...
    ,'parent',loadWellsPanel2,'fontsize',16,'fontunits','normalized','HorizontalAlignment','right'... % ,'HorizontalAlignment','left'
    ,'backgroundcolor','white','position',[0.2,0.86,0.3,0.06]);

% number of headers in file to load
uicontrol('parent',loadWellsPanel2,'style','edit','string',num2str(hDV.data.realWellData.loadedHeaderLines),...
    'pos',[.51 .87 .05 .06],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'callback',{@nWellHeaderscall});


% make uitable panel
uitableWellData = uitable('parent',loadWellsPanel2,'units','normalized','Position',[0.01 0.1 0.98 .72],'ColumnName', hDV.data.realWellData.inputStringColumnNames,...
    'data',hDV.data.realWellData.stringsWellDataAdvanced,'ColumnEditable', [false(size(hDV.data.realWellData.stringsWellDataAdvanced,2),1)]',...
    'tag','uitableWellDataTag1','CellEditCallback',{@editedUITableWellData});

%resize column name header font size
hs = '<html><h2>'; %html start
he = '</h2></html>'; %html end
cn = get(uitableWellData,'ColumnName'); %original
cnh = cellfun(@(x)[hs x he],cn,'uni',false); %with html
set(uitableWellData,'ColumnName',cnh) %apply



% %What's an API number button?
  pd=0.02; 
% uicontrol('parent',loadWellsPanel2,'style','pushbutton','string','API# help',...
%     'pos',[.2 pd .2 .1-2*pd],'fontsize',14,'fontunits','normalized',...
%     'callback',{@apiHelpButtonCallback});

%help format well data button
uicontrol('parent',loadWellsPanel2,'style','push','string','File Format Help',...
    'position',[0 pd .2 .1-2*pd],'fontsize',14,'fontunits','normalized','callback',{@callbuthelpformatdata,loadWellsFigure});

%text
uicontrol('Style','text','String',['Accepts up to ',num2str(hDV.data.nwells_max),' wells'] ... % ;' API#, X [km], Y [km], Year, month (1-12), Injection volume [bbl/month], Pressure [psi], Well name, and well type' }...
    ,'parent',loadWellsPanel2,'fontsize',16,'fontunits','normalized','HorizontalAlignment','right'... % ,'HorizontalAlignment','left'
    ,'backgroundcolor','white','position',[.65 pd .33 .1-2*pd]);

% checkbox: extrapolate injection to end of model or not?
if isfield(hDV.data.realWellData,'extrapolateInjectionCheck') % see if this variable exists, if it doesn't, make it. 
checkBoxValue=hDV.data.realWellData.extrapolateInjectionCheck;
else
    hDV.data.realWellData.extrapolateInjectionCheck=1;
   checkBoxValue=hDV.data.realWellData.extrapolateInjectionCheck ;
end
textwidth2=.25;
extrapInjectionCheckbox=uicontrol('parent',loadWellsPanel2,'style','checkbox','position',[.4+textwidth2 pd .05 .1-2*pd],...
                'backgroundcolor',[1 1 1],'callback',{@callbackExtapolateInjectionRate,hDV},'value',checkBoxValue);
 uicontrol('parent',loadWellsPanel2,'style','text','string','Extrapolate Injection?','position',[.4 pd-.01 textwidth2 .1-2*pd],...
                'fontsize',14,'fontunits','normalized','backgroundcolor',[1 1 1]);
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dropDownMenuWells([],[]); %start with the right panel on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%OK button
uicontrol('parent',loadWellsFigure,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',18,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,get(dropDownMenu2,'Value')});

valsout=0;

    function crf(~,~)
        hDV.tempvec=[];
        delete(loadWellsFigure);
        modal(hfig,'on');
    end

%this must tie back to callbackdatabuts.m
    function butcall(~,~,hDV,htxt,~)
        [P,Q]=size(htxt); vec=zeros(P,Q);
        val=get(dropDownMenu2,'Value'); % constant wells, or loaded well data?
        switch val
            case 1 % constant wells rate
                for p=1:P
                    for q=1:Q
                        vec(p,q)=str2double(get(htxt(p,q),'string'));
                    end
                end
                
                if  ~all(hDV.data.inject.vals==vec)
                   resetButtonsRed(hDV,[1,3,4])
                end
                hDV.data.inject.vals=vec;
                
                % make datenumBarrelsPerDay for constant rate wells.
%                    YearEndOfModel=get(hDV.hdsldr(1),'max');% end of scaleruler.
                % desired format for constant rate wells: thisWellDatenumBarrelsPerDay =
                % [ start_datenumber_of_injection,   0; ...
                %    end_datenumber_of_injection,   well_rate_bbls_per_day] ;
                for ffkj5=1:hDV.data.nwells
                    hDV.data.inject.datenumBarrelsPerDay{ffkj5,1}=[datenum(hDV.data.inject.vals(ffkj5,4),1,1,0,1,0), 0; datenum(hDV.data.inject.vals(ffkj5,5),1,1,0,1,0),hDV.data.inject.vals(ffkj5,3)] ;
                end

                hDV.data.inject.datenumBarrelsPerDay=hDV.data.inject.datenumBarrelsPerDay(1:hDV.data.nwells,1);    % clear any old wells out
                hDV.data.realWellData.use = 0 ;  %is it being used?
                
            case 2 % variable rate wells
                hDV.data.realWellData.stringsWellDataAdvanced = get(uitableWellData,'data');
                if ~isempty(hDV.data.realWellData.stringsWellDataAdvanced) % don't call function if no data is loaded in
                    SpreadsheetStrings2WellData(get(uitableWellData,'data'),hDV); % function that reads uitable data and saves it to structure
                else
                  msgWindow1=msgbox({cat(2,'You didn''t load any injection well data. If you don''t have any, manually enter at least one well, as you need one! ')}, 'Loading Well rate data','warn');
                  centerFigure(hDV.hfig,msgWindow1);
                  hDV.data.realWellData.use=0;
                end
                hDV.data.realWellData.use = 1 ;  %csv entry is being used
                hDV.data.inject.vals=zeros(P,Q); %this needs to be set somehow
        end

        delete(loadWellsFigure);
        modal(hfig,'on');
    end % end butcall function


    function nwellscall(src,~,hvalstxt,maxval)
        
        val = str2double(get(src,'string')) ;
        val=floor(val);
        if val<1
            val=1;
        elseif val>maxval
            val=maxval;
        end
        set(src,'string',num2str(val));
        
        hDV.data.nwells = str2double(get(src,'string')) ;
        
        for j=1:hDV.data.nwells
            set(hvalstxt(j,:),'visible','on') ;
        end
        
        for j=hDV.data.nwells+1:M
            set(hvalstxt(j,:),'visible','off') ;
        end
    end

    function dropDownMenuWells(~,~)  % based on dropdown menu, toggle visibility
        
        index_selected = get(dropDownMenu2,'Value');
        
        switch index_selected
            case 1 % Enter wells Manually
%                 disp('enter wells manually chosen')
                set(enterWellsPanel,'visible','on')
                set(loadWellsPanel2,'visible','off')
            case 2 % csv import
%                 disp('enter wells csv chosen')
                set(enterWellsPanel,'visible','off')
                set(loadWellsPanel2,'visible','on')
        end
    end


% % API help button
%     function apiHelpButtonCallback(~,~) % this callback opens the wikipedia page on API Well number in your default browser
%         url = 'https://en.wikipedia.org/wiki/API_well_number';
%         web(url,'-browser')
%     end


% help formatting well data image function
    function callbuthelpformatdata(~,~,hfig)
        modal(hfig,'off');
        hpos = get(hfig,'position') ;
        helpImageFigure = figure('color',[1 1 1],'units','pixels','visible','on');
        pos = [0 0 .6*hpos(3) .8*hpos(4)];
        imshow('wellDataLoadHelp.png','border','tight');
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

% callback reset headers
    function nWellHeaderscall(src,~)
        val = str2double(get(src,'string')) ;
        val=floor(val);
        if val<0
            val=0;
        end
        set(src,'string',num2str(val));
        hDV.data.realWellData.loadedHeaderLines=val;
    end



% read in csv file with monthly well data
    function callbackloadadvancedwellsfile(~,~) % read well data csv       
        [fileNameToLoad,pathname]=uigetfile('*.csv');
        fid44=fopen([pathname,fileNameToLoad]); % open file        
        readInData=textscan(fid44,'%s%s%s%s%s%s','Delimiter',',','Headerlines',hDV.data.realWellData.loadedHeaderLines); % read all data in as text strings
        fclose(fid44) ; % close file
        
        hDV.data.realWellData.stringsWellDataAdvanced=[readInData{:,:}];
        
        [hDV]=checkWellData(hDV);
        
        set(uitableWellData,'data',hDV.data.realWellData.stringsWellDataAdvanced)
        set(uitableWellData,'ColumnEditable',[false(size(hDV.data.realWellData.stringsWellDataAdvanced,2),1)]')
        
    end


% callback for editing uitable of well data
    function editedUITableWellData(hObject,callbackdata)
        %               editedVal = eval(callbackdata.EditData);
        editedVal = callbackdata.EditData;% new value
        oldVal=callbackdata.PreviousData; % old value
        r = callbackdata.Indices(1); % row edited
        c = callbackdata.Indices(2); % column edited
        [editedValasNumber, numericFlag] = str2num(editedVal); % try and convert to number, if unsuccessful, then numericFlag=0
        d=get(hObject,'Data');
        if numericFlag && hDV.data.realWellData.columnIsNumber(c)  % you entered numbers and column takes numbers
            % editedValasNumber = eval(callbackdata.EditData);
            d{r,c}=num2str(editedValasNumber);
            set(hObject,'Data',d);
        elseif ~numericFlag && ~hDV.data.realWellData.columnIsNumber(c) % you entered characters and column takes characters
            % then you entered characters
            d{r,c}=editedVal;
            set(hObject,'Data',d);
        else % error
            edlgBox=errordlg('you entered a character in a number column, or a number in a character column: edit rejected');centerFigure(hDV.hfig,edlgBox);
            % change data back to old data
            d{r,c}=oldVal;
            set(hObject,'data',d)
        end
    end

    function callbackExtapolateInjectionRate(~,~,~)
% callback function for if the "should we extrapolate a well's last injection rate to the end of the model?"
% check box is switched. 
% this is called at the end of SpreadsheetStrings2WellData
        hDV.data.realWellData.extrapolateInjectionCheck =get(extrapInjectionCheckbox,'Value');
        
       
    end

end




