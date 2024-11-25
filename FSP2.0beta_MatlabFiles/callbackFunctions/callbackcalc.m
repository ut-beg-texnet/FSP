% Callback when the calculate button is pushed

function callbackcalc(src,evt,hDV)

flag = checkdata(hDV) ; %check

if hDV.plotdata.printFunctionName
    disp(['running callbackcalc '])
end

%color of the calculate button
if flag
       % keep back compatible         If not designating seed
    if length(hDV.data.adv.vals)<=9 || hDV.data.adv.vals(10)==0
        rng('shuffle') % this resets the seed based on the computer clock each time
    else
        rng(round(hDV.data.adv.vals(10))) % set pseudorandom seed
    end
    %enable saving session
    set(hDV.bSave,'enable','on') ;
    %set the horizontal stress direction to correct mod
    hDV.data.stress.vals(4) = mod( hDV.data.stress.vals(4) , 180) ;
    
    if hDV.data.reservoir.importHydrology==0 % if using well data, not imported hydrologic model
        % decide if monthly well data, or constant well data
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
            minyrs=zeros(hDV.data.nwells,1);
            for k=1:1:hDV.data.nwells
                minyrs(k) = min(year(hDV.data.realWellData.datenumBarrelsPerDay{k}(:,1)));
                maxyrs(k) = max(year(hDV.data.realWellData.datenumBarrelsPerDay{k}(:,1)));
                
            end
            xmin = min(hDV.data.realWellData.XEasting);
            xmax = max(hDV.data.realWellData.XEasting);
            ymin = min(hDV.data.realWellData.YNorthing);
            ymax = max(hDV.data.realWellData.YNorthing);
            
        else % entered injection well constant rate
            %set the slider ranges
            minyrs=hDV.data.inject.vals(1:hDV.data.nwells,4); minyrs = minyrs(minyrs>0);
            maxyrs=hDV.data.inject.vals(1:hDV.data.nwells,5); maxyrs = maxyrs(maxyrs>0);
            %setup view grid auto
            xmin = min(hDV.data.inject.vals(1:hDV.data.nwells,1));
            xmax = max(hDV.data.inject.vals(1:hDV.data.nwells,1));
            ymin = min(hDV.data.inject.vals(1:hDV.data.nwells,2));
            ymax = max(hDV.data.inject.vals(1:hDV.data.nwells,2));
        end
        maxval=max([max(maxyrs);0])+2; %30 years past the max of well start year
        minval = min([min(minyrs);maxval-3 ]) ; 
        dl=maxval-minval;
        SliderStep=[1/dl 3/dl];
        grd=[xmin xmax ymin ymax] ;off=8;
        grd = grd+off*[-1 1 -1 1] ;
        
    else % hDV.data.reservoir.importHydrology ==1; % use imported hydrology
        yearsRepresentedHydroImport=hDV.data.reservoir.yearsRepresentedHydroImport;
        minval = min(yearsRepresentedHydroImport); maxval=max(yearsRepresentedHydroImport);
        if maxval==minval % if only 1 year of hydrology data entered:
            maxval=maxval+0.0001;
            SliderStep=[1,1];
        else % multiple years of hydrology imported
            SliderStep=[1/length(yearsRepresentedHydroImport),2/length(yearsRepresentedHydroImport)];
            
        end
    
        xmin=min(hDV.data.reservoir.numbersImportedHydrology(:,1));
        xmax=max(hDV.data.reservoir.numbersImportedHydrology(:,1));
        ymin=min(hDV.data.reservoir.numbersImportedHydrology(:,2));
        ymax=max(hDV.data.reservoir.numbersImportedHydrology(:,2));
        hDV.data.adv.vals(1:4)=[xmin xmax  ymin ymax];
        grd=[xmin xmax ymin ymax];
    end
        
    set(hDV.hdsldr(:),'min',minval,'max',maxval,'value',minval, 'SliderStep',SliderStep) ;
    set(hDV.hdsldr_txt(:),'string',num2str(minval)) ;

        
        
    switch hDV.data.fault.intype
        case 1
            %load faults parameter limits for random distributions
            nfaults = hDV.data.fault.vals(1); %number of faults
            Xmin=hDV.data.adv.vals(1); Xmax=hDV.data.adv.vals(2);
            Ymin=hDV.data.adv.vals(3); Ymax=hDV.data.adv.vals(4);
            Tmin=hDV.data.fault.vals(3);  Tmax=hDV.data.fault.vals(4);
            Dmin=hDV.data.fault.vals(5);  Dmax=hDV.data.fault.vals(6);
            mu=hDV.data.fault.vals(2);
                
            %randomize on each click
            buffer=(Xmax-Xmin).*0.05; % to make sure faults not at edge of map
            hDV.data.fault.xf=Xmin+buffer+(Xmax-Xmin-2*buffer)*rand(nfaults,1);
            hDV.data.fault.yf=Ymin+buffer+(Ymax-Ymin-2*buffer)*rand(nfaults,1);
            hDV.data.fault.thf=Tmin+(Tmax-Tmin)*rand(nfaults,1); %fault orientation
            hDV.data.fault.dipf =Dmin+(Dmax-Dmin)*rand(nfaults,1); %fault dips
            hDV.data.fault.lenf=4+zeros(nfaults,1) ;
            hDV.data.fault.muf=mu + zeros(nfaults,1) ;
                
                
        case 2    
            % columns are now: X easting, Ynorthing, strike, dip, length
            if iscell(hDV.data.fault.file)
                dat = cell2mat(hDV.data.fault.file);
            else
                dat = hDV.data.fault.file  ;
            end
            [nfaults,~]=size(dat) ;
            hDV.data.fault.xf=dat(:,1);
            hDV.data.fault.yf=dat(:,2);
            hDV.data.fault.thf=dat(:,3); %fault orientation deg
            hDV.data.fault.dipf =dat(:,4); %fault dips
            hDV.data.fault.lenf=dat(:,5);
            hDV.data.fault.muf=hDV.data.fault.vals(2)*ones(size(dat,1),1);  % All fault mu now have the same value FSP 2.0
            hDV.data.fault.vals(1) = nfaults ;
                
            if min(hDV.data.fault.xf)<grd(1), grd(1)=min(hDV.data.fault.xf)-1; end % this may need to be redone?
            if max(hDV.data.fault.xf)>grd(2), grd(2)=max(hDV.data.fault.xf)+1; end
            if min(hDV.data.fault.yf)<grd(3), grd(3)=min(hDV.data.fault.yf)-1; end
            if max(hDV.data.fault.yf)>grd(4), grd(4)=max(hDV.data.fault.yf)+1; end
                
    end
       
        
        %set the color bar range
        dpth=hDV.data.stress.vals(5); %calculation depth
        hDV.plotdata.mincfc=0 ; hDV.plotdata.maxcfc=round(max([min((hDV.data.stress.vals(1:3)-hDV.data.stress.vals(6)).*dpth)*1.2,200])); % 20% past frac gradient is Green end, or 200 PSI, whichever is larger
        set(hDV.plotdata.pffot.cmintxt,'string',num2str(hDV.plotdata.mincfc,'%5.2f'));
        set(hDV.plotdata.pffot.cmaxtxt,'string',num2str(hDV.plotdata.maxcfc));
        
        %reset the fault view
        set(hDV.plotdata.ListboxFaultSelector,'value',1)
%         set(hDV.plotdata.pffot.pop(2),'value',1); callfltpop(hDV.plotdata.pffot.pop(2),[],hDV) ;
        set(hDV.plotdata.inputMap.popWells,'value',1); callwellhighlight(hDV.plotdata.inputMap.popWells,[],hDV) ;
%         set(hDV.plotdata.pint.pop(2),'value',1); % callintfltpop(hDV.plotdata.pffot.pop(2),[],hDV) ;
        
        
        hDV.data.adv.vals(1:4)=grd; %set up data xmin/max ymin/max window
        resetprobpanel(hDV) ;
        callbackplots([],[],hDV,hDV.currtab.number,hDV.currtab.name);
        
        set(hDV.bCalc,'backgroundcolor', hDV.colors.green);
        set(hDV.bCalc,'foregroundcolor', hDV.colors.white);
        
else % flag
    
    set(hDV.bCalc,'backgroundcolor', hDV.colors.red);
    set(hDV.bCalc,'foregroundcolor', hDV.colors.white);
    %disable top buttons
    set(hDV.bMain.bTop(:),'Enable','off');
    
end % flag
    
        
end % function callbackcalc