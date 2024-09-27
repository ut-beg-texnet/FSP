% back door during setup to increase well or fault limit
% % include a .mat file in startup directory
% If the .mat file isn't in the same directory as the program, it won't be
% loaded
% running these 3 lines in matlab and putting this in the main directory
% edits the variables
% nwells_max=200; % 100 is default maximum
% NFAULTSMAX=600; % 500 faults is default maximum
% axsBgrnd = [1 , 1, 1] .*.4 % .8 is default RGB color
% mohrCircleDotSize = 4;  % 10 is default dot size
% save('./FSPStartup.mat','nwells_max','NFAULTSMAX','axsBgrnd','mohrCircleDotSize')


% 
function hDV = checkForStartupMatFile(hDV)

c=dir('.');
b=struct2cell(c);% is this filename in current directory?
if any(ismember(b(1,:),'FSPStartup.mat'))
% if mat file exists, load it
load('./FSPStartup.mat')
end
if exist('nwells_max','var')
    hDV.data.nwells_max=nwells_max;
    disp(['adding ',num2str(nwells_max),' wells limit at startup'])
end

if exist('NFAULTSMAX','var')
    hDV.data.NFAULTSMAX=NFAULTSMAX;
        disp(['adding ',num2str(NFAULTSMAX),' faults limit at startup'])
end

if exist('axsBgrnd','var')
    hDV.colors.axsBgrnd=axsBgrnd;
        disp(['setting axis background RGB to  [',num2str(axsBgrnd),'] at startup'])
end

if exist('mohrCircleDotSize','var')
    hDV.plotdata.mohrCircleDotSize=mohrCircleDotSize;
        disp(['setting mohr Circle Dot Size to ',num2str(mohrCircleDotSize),' at startup, 10 is default'])
end


end








