% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface


% This function checks the input data to ensure consistency (no negative
% values, values within range, etc.). Additional elseif statements can be
% added as needed.

function flag = checkdata(hDV)

%data checks
flag = 1;

%%% STRESS DATA TAB CHECKS

if hDV.data.stress.aphi.use==0 && (hDV.data.stress.vals(1)<=0 || hDV.data.stress.vals(2)<=0 || hDV.data.stress.vals(3)<=0)
    edlgBox=errordlg('Check stress gradients','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.stress.aphi.use==0 && hDV.data.stress.vals(2)>hDV.data.stress.vals(3) %check horizontal stresses
    edlgBox=errordlg('Check horizontal stress gradients','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.stress.aphi.use>0 && hDV.data.stress.vals(1)<=0 
    edlgBox=errordlg('Check vertical stress gradient','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.stress.aphi.use>0 && (hDV.data.stress.aphi.vals(1)<0 || hDV.data.stress.aphi.vals(1)>3)
    edlgBox=errordlg('Check A-Phi parameter value','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.stress.aphi.use==12 && hDV.data.stress.vals(2)<=0
    edlgBox=errordlg('Check minimum horizontal stress gradient','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.stress.vals(4)<0 || hDV.data.stress.vals(4)>=360
    edlgBox=errordlg('Check horizontal stress direction in range [0,360)','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.stress.vals(5)<=0
    edlgBox=errordlg('Check reference depth','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.stress.vals(6)<=0
    edlgBox=errordlg('Check initial reservoir pressure gradient','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.stress.vals(1)>1.5  % sV_grad > 1.5 psi/ft is not realistic
    edlgBox=errordlg('Check vertical stress gradient value','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.stress.vals(3)>8    % sHmax_grad>8 not realistic; no separate check for shmin as it is < sHmax 
    if hDV.data.stress.aphi.use == 0
        edlgBox=errordlg('Check horizontal stress gradient value(s)','Data Error') ; centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    elseif hDV.data.stress.aphi.use == 12
        edlgBox=errordlg('Check minimum horizontal stress gradient value','Data Error') ; centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    else
        edlgBox=errordlg('Calculated horizontal stresses too high; check vertical stress gradient and friction coefficient values','Data Error') ; centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    end
    
%%% HYDROLOGY DATA TAB CHECKS
elseif  hDV.data.reservoir.vals(5)<0 || hDV.data.reservoir.vals(5)>0.5
    edlgBox=errordlg('Check Poisson Ratio range [0,0.5]','Data Error') ;centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.reservoir.vals(1)<=0 && hDV.data.reservoir.importHydrology==0
    edlgBox=errordlg('Check aquifer thickness','Data Error') ;centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif (hDV.data.reservoir.vals(2)<1 ||   hDV.data.reservoir.vals(2)>=100) && hDV.data.reservoir.importHydrology==0
    edlgBox=errordlg('Check Porosity range [1% to 99%]','Data Error') ;centerFigure(hDV.hfig,edlgBox);
    flag=0;
elseif hDV.data.reservoir.vals(3)<=0 && hDV.data.reservoir.importHydrology==0
    edlgBox=errordlg('Check permeability','Data Error') ;centerFigure(hDV.hfig,edlgBox);
    flag=0;
    
    
    %%% FAULT DATA TAB CHECKS
elseif hDV.data.fault.vals(1)<1 || hDV.data.fault.vals(1)>hDV.data.NFAULTSMAX
    hDV.data.fault.vals(1) = floor(hDV.data.fault.vals(1)); %convert to integer
    edlgBox=errordlg(['Check number of faults range [1,',num2str(hDV.data.NFAULTSMAX),')'],'Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif ~isfield(hDV.data.fault,'file') || (isnan(hDV.data.fault.file(1)) &&  hDV.data.fault.intype==2 ) %no data when you've chosen enter faults
    edlgBox=errordlg('Check fault input','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.fault.vals(2)<=0 
    edlgBox=errordlg('Check friction coefficient for faults <=0','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.fault.vals(4)<hDV.data.fault.vals(3)
    edlgBox=errordlg('Check fault strike ranges','Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
elseif hDV.data.fault.vals(6)<hDV.data.fault.vals(5) || any(hDV.data.fault.vals(5:6)<0) || any(hDV.data.fault.vals(5:6)>90)
    edlgBox=errordlg(cat(2,'Check fault dip ranges: min ',num2str(hDV.data.fault.vals(5)),' and max ',num2str(hDV.data.fault.vals(6)) ),'Data Error') ; centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
    

%%% INJECTION WELLS TAB CHECK -- only if hydrology data is NOT imported
%%% from a reservoir model
elseif hDV.data.reservoir.importHydrology==0
    if ~isfield( hDV.data.realWellData,'use')
        edlgBox=errordlg('Check wells entry','Data Error') ;centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    elseif (~hDV.data.realWellData.use && all( hDV.data.inject.vals(:)==0))
        edlgBox=errordlg('Check Injection Wells entry','Data Error') ;centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    elseif (hDV.data.realWellData.use==1 && ~isfield(hDV.data.realWellData,'datenumBarrelsPerDay'))
         edlgBox=errordlg(['Check injection well entry, did you choose to load a .csv file but not enter any data? '],'Data Error') ;centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    elseif (~hDV.data.realWellData.use && all( hDV.data.inject.vals(1:hDV.data.nwells,3)==0))
        edlgBox=errordlg('Check Injection Wells entry: did you enter 0 injection rate wells?','Data Error') ;centerFigure(hDV.hfig,edlgBox);
        flag=0 ;
    end

%%% ADVANCED Tab Checks
elseif hDV.data.adv.vals(10)~=round(hDV.data.adv.vals(10)) || hDV.data.adv.vals(10)<0
    edlgBox=errordlg(['Check Seed Set in Advanced tab=',num2str(hDV.data.adv.vals(10)),' , 0 doesn''t set a seed, a positive integer sets the seed'],'Data Error') ;centerFigure(hDV.hfig,edlgBox);
    flag=0 ;
    
    % warn if large area? 
end






