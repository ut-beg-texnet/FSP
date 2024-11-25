function callbackloadsess(~,~,hDV)

uiopen('LOAD'); 

if exist('chkstring','var')
if strcmpi(chkstring,'DARREN PAIS URC 2016') || strcmpi(chkstring,'Saved FSP Session');
    RunningVersion=hDV.data.version;
    disp('File Loaded') ;
    hDV.data =datsviz;
    num=1 ; callbackplots([],[],hDV,num,hDV.tabnames{num}) ;
    set(hDV.bCalc,'enable','on') ;
    
    % backwards compatible from 0.94 
    if ~isfield(hDV.data.reservoir,'importHydrology')
        hDV.data.reservoir.importHydrology=1;
    end
    if ~isfield(hDV.data.reservoir,'loadedHeaderLines')
        hDV.data.reservoir.loadedHeaderLines=1;
    end
    if ~strcmpi(hDV.data.version,RunningVersion)
    % heads up if different versions
    msgWindow1=msgbox(['You just loaded an FSP file that was saved with version ',hDV.data.version,' but you are running version ',...
        RunningVersion, ', the session will load, however some plots may appear different'], 'FSP version warning','warn');
    centerFigure(hDV.hfig,msgWindow1);
    end
    if ~isfield(hDV.data,'NFAULTSMAX') % 0.99.2 moved this. Should be able to get rid of after 1.0
        %         hDV.data.NFAULTSMAX=hDV.NFAULTSMAX;
        hDV.data.nwells_max = hDV.plotdata.wellsMaxFaultsMax(1) ;
    end
    
    if ~isfield(hDV.data,'nwells_max') % 0.99.2 moved this. Should be able to get rid of after 1.0
        hDV.data.nwells_max = hDV.plotdata.wellsMaxFaultsMax(1) ;
    end

    if length(hDV.data.adv.vals)==9 % for before added 
        hDV.data.adv.vals(10)=0;
        disp(['setting hDV.data.adv.vals(10)=0 because not in loaded data'])
    end
    if ~isfield(hDV.data.reservoir,'importHydrology') % 0.99.2 moved this. Should be able to get rid of after 1.0
      disp(['setting hDV.data.reservoir.importHydrology=0 because not in loaded data'])
        hDV.data.reservoir.importHydrology=0;
    end
    if ~isfield(hDV.data.realWellData,'use') ;
    hDV.data.realWellData.use =  0 ;
    end
else
    msgWindow1=msgbox(['Check String Load Error, are you sure this is an FSP saved session?'], 'FSP String Load Error','warn');
    centerFigure(hDV.hfig,msgWindow1);
end
else
    disp('No File Loaded') ; 
     msgWindow1=msgbox(['No File was loaded'], 'FSP file Load Error','error');
    centerFigure(hDV.hfig,msgWindow1);
    
end
end