classdef surfaceviz < DASVdispQRA
    
    %--------------------------------------------------------------------------
    % (c) ExxonMobil Upstream Research Company, 2014, 2016++
    %     Drilling & Subsurface Technology Division
    %     Well Construction Section
    %  Authors: Payette, GS; Edited: Pais, D.
    %  GUI Setup for QRA tool
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    %
    % Class properties
    %
    %--------------------------------------------------------------------------
    
    %----------------------------Public properties----------------------------%
    properties (Access = 'public')
        
        % Primary panels and buttons + data structures
        pCntr      = [];    % controls panel handle
        pMain      = [];    % main panel handles
        hExit      = [];    % handle to exit button
        bCntr      = [];    % controls panel buttons and data
        bMain      = [];    % data structure with main panel(s) buttons
        bCalc      = [];    % calculate button
        bDef       = [];    % default data button
        bSave      = [];    % save session button
        bLegal     = [];    % legal disclaimer button
        bLoad      = [];    % load session button
        tabnames   = [];    % "tab" names for each analysis module
%         cbutnames  = [];    % names for control button % now under uimenuhandles
        currtab    = [];    % current selected calculation "tab"
        hdsldr     = [];    % handle to the time slider object (this is a vector of all identical time sliders)
        hdsldr_txt = [];    % vector of all idential time slider text handles
        colors     = [];    % colors; used for buttons
        ftsz       = 14;    % fontsize for axes
        cmapGYR    = [];    % Green Yellow Red colormap
        cmapInfo   = [];    % For horizontal color bars
        tempvec    = [];    %temporary storage vector
        uiMenuHandles=[];   % menu handles... File, Zoom, data inputs, etc.  
        
        % Directory, software version, reference figure / button
        hfigMain   = [];    % handle to main figure
        hbutMain   = [];    % handle to button on main figure
        
        
        % Data handle
        data       = [];    % input data goes here
        plotdata   = [];    % data arrays for plots go here
        
%         NFAULTSMAX = 500;
    end
    
    %--------------------------------------------------------------------------
    %
    % Public methods: class constructor, destructor and accessor functions
    %
    %--------------------------------------------------------------------------
    
    methods (Access = 'public')
        
        %----------------------------Class constructor----------------------------%
        function hDV = surfaceviz(vers)
            
            %-------------------------Preliminary routines-------------------------%
            
            % Call the DASdisp constructor to make the window
            hDV = hDV@DASVdispQRA(['Fault Slip Potential v',vers]);
            
            
            
            % Colors
            hDV.colors.white     = [1 1 1];
            hDV.colors.black     = [0 0 0];
            hDV.colors.red       = [1 0 0];
            hDV.colors.darkred   = [0.8 0 0];
            hDV.colors.grey      = 0.8 * [1 1 1];
            hDV.colors.darkgrey  = 0.3 * [1 1 1];
            hDV.colors.default   = 0.941176 * [1 1 1]; % Matlab "special" default
            hDV.colors.blue      = [0 0 0.8];
            hDV.colors.yellow    = [1 1 0];
            hDV.colors.green     = [0 .7 0];
            hDV.colors.axsBgrnd  = [1 1 1] .*.8; % axis background so you see yellows
            
            %colormap stuff
            x= load('cmaps','cmapGYR');
            hDV.cmapGYR = x.cmapGYR ;
            hDV.cmapInfo.patchesInColorbar=300;
            
            %max number of wells
            hDV.data.nwells_max = 100 ; 
            hDV.data.NFAULTSMAX = 500;
            hDV.data.version=vers;
            hDV.plotdata.mohrCircleDotSize = 10;
            hDV = checkForStartupTxtFile(hDV); % back door to increase fault or well count limit
            hDV.plotdata.wellsMaxFaultsMax=[hDV.data.nwells_max;hDV.data.NFAULTSMAX];
            %----------Setup Fault Slip Potential logo image, primary panels and buttons----------%
            
            % Setup the "FSP logo" image
            ax1 = axes('parent',hDV.hfig,'position',[0.015 0.895 0.15 0.13]);
            rgb = imread('logo.png');
            image(rgb,'parent',ax1);
            axis(ax1,'off','image');
            
            setupdata(hDV) ; %setup the data array
            
            
            % Parameters for setting panel locations on main fig
            GW = 0.005;             % gutter width
            HP = 0.86;              % height of the two panels
            H1 = 0.07;              % distance from fig bottom to panel bottom
            WC = 0.18;              % width
            WM = 1 - (WC + 3 * GW); % width of main panel
            HE = 0.06;              % exit button height
            WE = 0.15;              % exit button width
            
            % Positions: controls panel, main panel and exit button
            posCntr = [GW H1 WC HP];
            posMain = [(1 - GW - WM) 0.5*(H1 - HE) WM HP+.065];
%             posExit = [(GW + 0.5*(WC - WE)) 0.5 * (H1 - HE) WE HE];
            
            % Controls panels
            hDV.pCntr = uipanel('parent',hDV.hfig,'backgroundcolor','white',    ...
                'visible','on','Position',posCntr);
            
            %calculate button
            hDV.bCalc = uicontrol('parent',hDV.pCntr,'Style','push','String','Calculate',...
                'pos',[.05 .02 .9 .05],'fontsize',12,'fontunits','normalized',...
                'TooltipString',sprintf('%s','Run calculations after you''ve entered in your data'),...
                'callback',{@callbackcalc,hDV},'enable','on');
            
           
            bottom=get(hDV.bCalc,'pos');
            
            uicontrol('parent',hDV.pCntr,'style','text','string','Fault Selector:',...
                'pos',[.1 .95  .8 .05],'fontsize',12,'fontunits','normalized',...
                'backgroundcolor',[1 1 1]);
            
            top=.95;
                       
            % multiple fault selector menu
            space=0.02;
            bottom=bottom(2)+bottom(4)+space;
            
            positionListbox=[.05 bottom .9 top-bottom];
            hDV.plotdata.ListboxFaultSelector=uicontrol('style','listbox','parent',hDV.pCntr,'string',{'Select faults here'},'position',positionListbox,...
                'fontsize',12,'fontunits','normalized','value',1,'backgroundcolor',[1 1 1],...
                'min',1,'max',10,'callback',{@callbackFaultsSelected,hDV},...
                'TooltipString','hold control to select multiple faults, hold shift to select multiple faults in a row');
           
            %tabs to use
            hDV.tabnames = {'MODEL INPUTS' , 'GEOMECHANICS' , 'PROB. GEOMECH' , 'HYDROLOGY' ,'PROB. HYDRO' , 'INTEGRATED'};
            num       = length(hDV.tabnames);
            % yellow button explanation when moused over but not clicked
            tabTooltipStrings={'Map of Faults and Wells, but no analysis';'Deterministic Geomechanical Analysis';'Probabilistic Geomechanical Analysis of Pore Pressure to Slip';...
                'Deterministic Hydrological Analysis';'Probabilistic Hydrological analysis of Pressure on Faults';'Integrating Fault Slip Potential'};
            
            % Parameters for setting button locations on top of main fig
            WB  = (WM - GW * (num + 1)) / num;      % width of each button
            HB  = (1 - (H1 + HP + 5 * GW));         % height of each button
            pos = [(3 * GW + WC) (H1 + HP) WB HB];  % position of first button
            del = [(WB + GW) 0 0 0];                % translation of other buttons
            
            % Heights of buttons at top of screen (non-selected and selected)
            hDV.bMain.bTopHeights.HB    = HB;
            hDV.bMain.bTopHeights.HBdel = 4*GW;
            
            % Setup the tab buttons at the top
            hDV.bMain.bTop = NaN * zeros(num,1);
            for i = 1:1:num
                hDV.bMain.bTop(i) = uicontrol('parent',hDV.hfig,'Style','push',  ...
                    'String',hDV.tabnames{i},'pos',pos + (i-1)*del,               ...
                    'HandleVisibility','off','fontsize',12,'fontunits',            ...
                    'normalized','callback',{@callbackplots,hDV,i,hDV.tabnames{i}}) ;
                set(hDV.bMain.bTop(i),'TooltipString',sprintf('%s',tabTooltipStrings{i})); % yellow button explanation when moused over but not clicked
            end
            
            % Main panel(s): multiple panels defined for easy linking with tabs
            hDV.pMain = NaN * zeros(num+1,1);   % make one extra "blank" panel
            for i = 1:1:(num+1)
                hDV.pMain(i) = uipanel('parent',hDV.hfig,'backgroundcolor',      ...
                    'white','visible','off','Position',posMain);
            end
            set(hDV.pMain(end),'visible','on');
            setupplotpanels(hDV) ; % setup all the plotting panels
            
            % make uimenu dropdowns 
            makeFileUiMenu(hDV) % file, quit, save, load, etc menu
            makeInputDataUiMenu(hDV)
            makeExportUImenu(hDV) % menu of which axes image can be exported
            makeZoomUiMenu(hDV)
            
            %default start panel
            hDV.currtab.name=hDV.tabnames{1} ;
            hDV.currtab.number =1 ;
            
            %disable top buttons
            set(hDV.bMain.bTop(:),'Enable','off');
            
            % Setup the buttons and sliders for time windowing
            %%% post SCITS edit: Date Slider only embeded in selected panels
            %PI = [(1 - WM - GW) 0.005];
            %dateSldr(hDV.hfig,PI,WM,HE,GW,hDV); %create date slider
            
        end  % surfaceviz method (class constructor)
        
        %----------------------------Class destructor-----------------------------%
        function delete(hDV)
            
            %put destructor stuff here
            
        end  % delete method (class destructor)
        
    end  % public methods
    
    %--------------------------------------------------------------------------
    %
    % Static methods: these methods may be run without ever creating an object
    %
    %--------------------------------------------------------------------------
    
    methods (Static)
        
        %-------------------------Write subtext in a box--------------------------%
        function hout = subtext(hparent,pos,text,fontsz,colr,colr2)
            if ( nargin == 5 )
                colr2 = 'cyan';
            end
            hout = uicontrol('parent',hparent,'Style','Text','Position',pos,    ...
                'String',text,'fontunits','normalized','backgroundcolor',colr,      ...
                'foregroundcolor',colr2,'fontsize',fontsz,'HorizontalAlignment',    ...
                'center');
        end
        
    end  % static methods
    
end % classdef

%--------------------------------------------------------------------------
%
% End of class: class dismissed
%
%--------------------------------------------------------------------------