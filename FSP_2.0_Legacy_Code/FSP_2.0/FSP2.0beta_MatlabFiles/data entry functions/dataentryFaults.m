% Stanford 2016
% Rall Walsh with code from Darren Pais
% possible additional feature: add blank row button next to help button.
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
pos = [0 0 .5*hpos(3) .8*hpos(4)];
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

radioButtonGroup=uibuttongroup('Position', [0.1 .75 .4 .1],'backgroundcolor',[1 1 1],...
    'SelectionChangeFcn', {@radioMenuFaultsToggle},'parent',f,'fontsize',hDV.ftsz,...
    'BorderType','none','fontunits','normalized');

radioButton1=uicontrol(radioButtonGroup,'Style','radiobutton',...
    'Position', [0.01 .45 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Random Faults','fontunits','normalized',...
    'ToolTipString','Select to randomly generate faults');

radioButton2=uicontrol(radioButtonGroup,'Style','radiobutton',...
    'Position', [0.01 .0 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Enter Faults','fontunits','normalized',...
    'ToolTipString','Select to type in faults or load a .csv file');


valtxtH = zeros(length(vals),1) ; % preallocate

% Create common boxes, number of faults and fric coeff, for both random and
% manually entered faults

uicontrol('parent',f,'style','text','string',txt{1},...
    'pos',[.05 .92  .45 .04],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);
        
valtxtH(1)= uicontrol('parent',f,'style','edit','string',num2str(vals(1)),...
    'pos',[.45 .92  .45 .04],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],...
    'callback',@callbackFaultCountField);

uicontrol('parent',f,'style','text','string',txt{2},...
    'pos',[.05 .86  .45 .04],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);
        
valtxtH(2)= uicontrol('parent',f,'style','edit','string',num2str(vals(2)),...
    'pos',[.45 .86  .45 .04],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);



% type faults table
enterFaultsPannel = uipanel('parent',f,'backgroundcolor','white',    ...
    'visible','off','Position',[0 0.1 1 .64]);

% uicontrol('parent',enterFaultsPannel,'style','text','string',txt{1},...
%     'pos', [.05 .94 .45 .05],'fontsize',12,'fontunits','normalized',...
%     'backgroundcolor',[1 1 1]);
faultTable = uitable('Parent', enterFaultsPannel,'units','Normalized' ,'fontunits','normalized','fontsize',.04, 'Position', [0.1 .1 .8 .9]);

% Controls panels
% make random faults pannel
randomFaultsPannel = uipanel('parent',f,'backgroundcolor','white',    ...
    'visible','on','bordertype','none','Position',[0 0.25 1 .53]);


% make random faults pannel
pos = [0.05 0 .45 .06] ;  c2=[.45 0 -.15 0] ;
m=1 ;
for k=3:1:length(txt)
    if ~any(k==noshow)
        uicontrol('parent',randomFaultsPannel,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.11 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        
        valtxtH(k)= uicontrol('parent',randomFaultsPannel,'style','edit','string',num2str(vals(k)),...
            'pos',c2+pos+[0 .95-m*.11 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        m=m+1 ;
    else
        valtxtH(k)= uicontrol('parent',randomFaultsPannel,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
end


% Set fault table length according to input of faults input
%
callbackFaultCountField(valtxtH(1),[])
cellRowsToPopulate = str2double(get(valtxtH(1),'string'));

% get data for table from structure if it exists
if ~isfield(hDV.data.fault,'file')
    indat = cell(cellRowsToPopulate,5);
else
    indat = cell(cellRowsToPopulate,5);
    
    if isnumeric(hDV.data.fault.file) && all(all(~isnan(hDV.data.fault.file)))
        vals = num2cell(hDV.data.fault.file);
        [M,N]=size(vals);
        if N>5
            N=5;  % Use only first 5 columns, ignore 6th for fric coeffs
        end
        indat(1:M,1:N)=vals(1:M,1:N);
        set(valtxtH(1),'string',num2str(M)); % set number of faults to match
        % if spreadsheet data was previously loaded, show that
        if hDV.data.fault.intype==2
            set(radioButtonGroup,'SelectedObject',radioButton2)
            set(randomFaultsPannel,'visible','off')
            set(enterFaultsPannel,'visible','on')
        end
    end
end


%  X[easting km], Y[northing km], strike [deg], dip [deg],length [km]
set(faultTable, 'Data', indat);
set(faultTable, 'ColumnName' , {'<html><h3>  &nbsp X [East km] &nbsp   </h3></html>', '<html><h3> &nbsp Y [North km] &nbsp</h3></html>', '<html><h3> &nbsp Strike [Deg]  &nbsp </h3></html>', '<html><h3>  &nbsp Dip [Deg] &nbsp</h3></html>',...
    '<html><h3> &nbsp Length [km] &nbsp </h3></html>'});
set(faultTable, 'ColumnWidth', {'auto', 'auto', 'auto', 'auto', 'auto' });
set(faultTable, 'ColumnEditable' , [true true true true true ]);


%load faults file button
uicontrol('parent',enterFaultsPannel,'Style','push','String','Load File',...
    'pos',[0.1,0,.3,.1],'fontsize',12,'fontunits','normalized',...
    'callback',{@callbackloadfaultsfile,faultTable,hfig,hDV});

%help  button
uicontrol('parent',enterFaultsPannel,'Style','push','String','Help',...
    'pos',[.4,0,.3,.1],'fontsize',12,'fontunits','normalized',...
    'callback',{@callbutHelpFormatFaultData,f});


%OK button
uicontrol('parent',f,'style','pushbutton','string','OK',...
    'pos',[.3 .02 .4 .06],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,radioButtonGroup});


% Return value for dataentryFaults
%
valsout=0;



    function crf(~,~)
        hDV.tempvec=[];
        delete(gcf);
        modal(hfig,'on');
    end


    % Callback function for "Load File" button to read in csv file
    %
    function callbackloadfaultsfile(~,~,t,hfig,hDV)
        [fileNameToLoad,pathname]=uigetfile('*.csv');
        nHeaderlinesCSVFault=0;
        readInData=[];
        NumColumnsDesired=5; % desired data column number
        minColumns=3; % minimum number of csv columns accepted
        
        
        assumedFaultLength=2; % assumend fault length (km), if not loaded
        assumedDip=80; % degrees
                
        
        if fileNameToLoad~=0
            try
                readInData=csvread([pathname,fileNameToLoad],nHeaderlinesCSVFault,0); % read in csv file of numbers
                [rowsReadIn,colsReadIn]=size(readInData);
                if colsReadIn>NumColumnsDesired % check number of columns
                    errorWindow1 = errordlg(cat(2,'you entered in ',num2str(colsReadIn),' columns of data in the .csv file, but  ',num2str(NumColumnsDesired),' are needed, throwing out any extra column(s)'));
                    centerFigure(hfig,errorWindow1);
                    readInData=readInData(:,1:NumColumnsDesired);
                elseif colsReadIn<minColumns % check number of columns
                    errorWindow1 = errordlg(cat(2,'you entered in ',num2str(colsReadIn),' columns of data in the .csv file, but netween ',num2str(minColumns),' and ',num2str(NumColumnsDesired),' are needed,'));
                    centerFigure(hfig,errorWindow1);
                elseif colsReadIn==4
                    errorWindow1 = errordlg(cat(2,'you didn''t enter fault length for any fault, so assuming each fault is ', num2str(assumedFaultLength), ' km long '));
                    centerFigure(hfig,errorWindow1);
                    readInData(:,5)=assumedFaultLength;
                elseif colsReadIn==3
                    errorWindow1 = errordlg(cat(2,'you didn''t enter fault length, or dip for any fault, so assuming each fault is ',...
                        num2str(assumedFaultLength), 'km long, and has a dip of  ',num2str(assumedDip),' degrees'));
                    centerFigure(hfig,errorWindow1);
                    readInData(:,4)=assumedDip;
                    readInData(:,5)=assumedFaultLength;
                end
                if rowsReadIn>hDV.data.NFAULTSMAX % check number of rows
                    errorWindow1= errordlg(cat(2,'you entered in ',num2str(rowsReadIn),' rows of data in the .csv file, but only ',num2str(hDV.data.NFAULTSMAX),...
                        ' are accepted, throwing out extra row(s) of data'));
                    centerFigure(hfig,errorWindow1);
                    readInData=readInData(1:hDV.data.NFAULTSMAX,:);
                end
                if any(any(isnan(readInData)))
                    [rowsNaN,colsNaN]=find(isnan(readInData));
                    allOneString = sprintf('%.0f,' , rowsNaN);
                    allOneString2 = sprintf('%.0f,' , colsNaN);
                    errorWindow1= errordlg(cat(2,'you entered something that is Not a Number in row(s) ',allOneString,' and column(s) ',allOneString2,...
                        ' (respectively). Fix this manually or re-load the csv file'));
                    centerFigure(hfig,errorWindow1);
                end
               
                
                set(t,'Data',readInData) % move data to uitable
                set(valtxtH(1),'string',num2str(rowsReadIn)) % update number of faults
                
            catch err
                if strcmpi('MATLAB:textscan:handleErrorAndShowInfo',err.identifier) % if anything entered is not a number
                    errorWindow1= errordlg({'you may have a letter in this csv file. ';err.message});
                    centerFigure(hfig,errorWindow1);
                else
                    rethrow(err);
                end
            end
        end
        
    end



    % Callback function to toggle visibility of panels based on radio button
    % 
    function radioMenuFaultsToggle(src,event)  
        
        index_selected = get(event.NewValue,'String'); % get(dropDownMenu1,'Value');
        
        switch index_selected
            case 'Random Faults' % Random Faults
                set(enterFaultsPannel,'visible','off')
                set(randomFaultsPannel,'visible','on')
            case 'Enter Faults' % Enter Faults
                set(randomFaultsPannel,'visible','off')
                set(enterFaultsPannel,'visible','on')
        end
        
    end



    %this must tie back to callbackdatabuts.m  (SPL: is this still relevant?)
    % Callback function for OK button
    %
    function butcall(~,~,hDV,htxt,radioButtonGroup)
        
        stringSelected= get(get(radioButtonGroup,'SelectedObject'),'string');
        switch stringSelected
            case 'Random Faults'
                inputMode=1;
            case  'Enter Faults'
                inputMode=2;
        end
        
        hDV.data.fault.intype = inputMode;
        
        vec=zeros(length(htxt),1);
        for n=1:1:length(htxt)
            vec(n)=str2double(get(htxt(n),'string'));
        end
        
        if inputMode==1  % random faults
            
            if ~all(hDV.data.fault.vals==vec)
                resetButtonsRed(hDV,[1,2,3,4]) % Set all buttons red
            end
                
            hDV.data.fault.vals=vec;
            hDV.data.fault.file= NaN;
                
            delete(gcf);
            modal(hfig,'on');
                
        elseif inputMode==2 % entered manually
            
            % Set number of faults and friction coefficients from text
            % boxes above the table in the input window 
            %
            hDV.data.fault.vals(1) = vec(1);
            hDV.data.fault.vals(2) = vec(2);
            
            hDV.data.fault.file= get(faultTable,'data');
            [~,cols]=size(hDV.data.fault.file);
            
            if iscell(hDV.data.fault.file) % If fault data entered manually in table
                
                if all(all(cellfun('isempty',hDV.data.fault.file)))
                    ed=errordlg('you didn''t enter any faults! switch to randomly generated faults, or enter fault data','Faults Entry','modal') ;
                    centerFigure(f,ed);
                    waitfor(ed);
                    hDV.data.fault.file = [];
                else % try to convert to matrix
                    try
                        hDV.data.fault.file= cell2mat(hDV.data.fault.file);
                    catch err
                        if strcmpi('MATLAB:catenate:dimensionMismatch',err.identifier) % if anything entered is not a number
                            errorWindow1= errordlg({'you may have incomplete rows of data entered ';err.message});
                            centerFigure(hfig,errorWindow1);
                        else
                            rethrow(err);
                        end
                        hDV.data.fault.file = [];
                    end    
                end
                
            end
            
            
            if  cols==1 && isnan(hDV.data.fault.file);
                ed=errordlg('you didn''t enter any faults! switch to randomly generated faults, or enter fault data','Faults Entry','modal') ;
                centerFigure(f,ed);
                waitfor(ed);
            elseif  cols~=5
                ed=errordlg('Check Data Entry for Faults, number of columns does not = 5','Faults Entry','modal') ;
                centerFigure(f,ed);
                waitfor(ed);
            elseif any(any(isnan(hDV.data.fault.file)))
                [rowsNaN,colsNaN]=find(isnan(hDV.data.fault.file));
                allOneString = sprintf('%.0f,' , rowsNaN);
                allOneString2 = sprintf('%.0f,' , colsNaN);
                errorWindow1= errordlg(cat(2,'you entered something that is Not a Number in row(s) ',allOneString,' and column(s) ',allOneString2,...
                    ' (respectively). Fix this manually or re-load the csv file before proceeding.'));
                centerFigure(hfig,errorWindow1);
            else               
                
                delete(f);
                modal(hfig,'on');
                resetButtonsRed(hDV)
            end
            
            % Set currently selected faults to first fault, if any have been
            % selected earlier
            %
            if isfield(hDV.plotdata, 'curfault')
                    hDV.plotdata.curfault = 0;
            end
             
        end % if inputmode
        
        
        % If APhi model is being used, update the corresponding horizontal
        % stress (only so that boxes on stress dataentry tab are pre-populated)
        % based on new friction coeff value
        % FSP 2.0 uses fault_mu as the mu for APhi model, so need to update these here
        %
        if hDV.data.stress.aphi.use==11  % i.e. only APhi, NOT Modified-APhi
            APhi = hDV.data.stress.aphi.vals(1);
            Sv   = hDV.data.stress.vals(1);
            pp   = hDV.data.stress.vals(6);
            aphi_use = hDV.data.stress.aphi.use;
            [SH,Sh] = getHorFromAPhi(APhi,vec(2),0.000,Sv,pp,aphi_use); % input #3, shmin is irrelevant here
            hDV.data.stress.vals(2) = Sh;
            hDV.data.stress.vals(3) = SH;
        end
               
    end % callback for ok button


    % Callback function for "Help" button to show formatting fault input data help image
    %
    function callbutHelpFormatFaultData(~,~,f)
        modal(f,'off');
        hpos = get(f,'position') ;
        f2 = figure('color',[1 1 1],'units','pixels','visible','on');
        pos = [0 0 .8*hpos(3) 1*hpos(4)];
        imshow('images/faultDataLoadHelp.png','border','tight');
        set(f2,'position',pos) ;
        set(f2,'MenuBar','none');
        set(f2,'NumberTitle','off');
        set(f2,'DefaultUicontrolUnits','normalized');
        set(f2,'DefaultUicontrolFontsize',14);
        set(f2,'closerequestfcn',{@crf,f,f2}) ;
        centerFigure(f,f2);
        %set(f,'visible','on') ;
        function crf(~,~,f,f2)
            delete(f2);
            modal(f,'on');
        end
        
    end


    % Callback function to adjust length of faulttable according to number of
    % faults input provided
    %
    function callbackFaultCountField(src,~) 
        
        currentTableCellData=get(faultTable,'data'); % current table data
        numCols=size(currentTableCellData,2); % number of columns in table (5)
        cellRowsToPopulate = str2double(get(src,'string')); % number of rows desired in table
        if isnan(cellRowsToPopulate) || cellRowsToPopulate<=0 % make sure bounds reasonable
            cellRowsToPopulate=20;
        elseif cellRowsToPopulate>hDV.data.NFAULTSMAX % if too high
            cellRowsToPopulate=hDV.data.NFAULTSMAX;
        else % round down
            cellRowsToPopulate=floor(cellRowsToPopulate);
        end
        set(src,'string',num2str(cellRowsToPopulate)) % set field to value
        % now change UItable size:
        if cellRowsToPopulate > size(currentTableCellData,1) % if adding rows to uitable
            newCells=cell(cellRowsToPopulate-size(currentTableCellData,1),numCols); % make blank cells
            currentTableCellData =  [currentTableCellData;newCells];
        elseif cellRowsToPopulate < size(currentTableCellData,1) % if subtracting rows from uitable
            currentTableCellData= currentTableCellData(1:cellRowsToPopulate,:);
        end
                
        set(faultTable,'data',currentTableCellData ) % make table size
        
    end
end