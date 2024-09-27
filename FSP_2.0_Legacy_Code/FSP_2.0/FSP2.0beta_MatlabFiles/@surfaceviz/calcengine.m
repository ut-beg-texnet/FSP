%calculation engine calls. This function sets up and runs the mohrs_3D.m
%
% Modified by Suvrat Lele for FSP 2.0
%

function calcengine(hDV,name)

if hDV.plotdata.printFunctionName
    disp(['running calcengine for ',name])
end
switch name
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'MODEL INPUTS'
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case {'HYDROLOGY','INTEGRATED'}
        
        if hDV.data.reservoir.importHydrology==0; % 0 does hydrology calculation internally, 1 imports XYZPT model like from MODFLOW
        % waitbar
            waitBarHandleWells = waitbar(0,'Calculating ...','Name','Hydrology Calculation');
            centerFigure(hDV.hfig,waitBarHandleWells)
        
        %%% field surface
        xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
        ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
        npts=50 ; nwells=hDV.data.nwells;
        [hDV.plotdata.pflot.Xgrid, hDV.plotdata.pflot.Ygrid] = meshgrid(linspace(xmin,xmax,npts) , linspace(ymin,ymax,npts)) ;
        r=linspace(0.5,20,npts); %km for pressure curves
        
         waitbar(.02,waitBarHandleWells,['Calculating ',num2str(2,'%3.0f'), '% ...']) ;
        hDV.plotdata.pflot.Zgrid = pfieldcalc(hDV); %in psi
        
        waitbar(1 / 3,waitBarHandleWells,['Calculating ',num2str(100*1 / 3,'%3.0f'), '% ...']) ;
        
        %%% pressure curves
        hDV.plotdata.pflot.x=zeros(nwells,length(r));
        hDV.plotdata.pflot.y=hDV.plotdata.pflot.x; % reinitialize data arrays
        
        for k=1:hDV.data.nwells
            %location of well
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                xwell=hDV.data.realWellData.XEasting(k);
                ywell=hDV.data.realWellData.YNorthing(k);
            else
                xwell=hDV.data.inject.vals(k,1); ywell=hDV.data.inject.vals(k,2);
            end
            
            p = pfieldcalc(hDV,xwell+r,ywell,k); %radial flow model assumed
            %update plot data
            hDV.plotdata.pflot.x(k,:)=r; hDV.plotdata.pflot.y(k,:)=p;
            fraction=(1/3)+(k/(3.*hDV.data.nwells));
            waitbar(fraction,waitBarHandleWells,['Calculating ',num2str(100*fraction,'%3.0f'), '% ...']) ;
        end
       
        
        hDV.plotdata.pint.Xgrid = hDV.plotdata.pflot.Xgrid;
        hDV.plotdata.pint.Ygrid = hDV.plotdata.pflot.Ygrid;
        hDV.plotdata.pint.Zgrid = hDV.plotdata.pflot.Zgrid;
        
        
        %fault failure calculation, Mohrs circle
        thf=hDV.data.fault.thf ;   %fault strikes
        dips=hDV.data.fault.dipf; %fault dips
        SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
        muf=hDV.data.fault.muf; %fault mu
        
        %get pressure from hydrology model at each fault
        %ppf = getP_interp(hDV.data.fault.xf,hDV.data.fault.yf,hDV.plotdata.pflot.Xgrid, hDV.plotdata.pflot.Ygrid,hDV.plotdata.pflot.Zgrid) ; %in psi
        nfaults = hDV.data.fault.vals(1) ; ppOnFault=zeros(nfaults,1) ;
        for j=1:1:nfaults
            wells=1:1:hDV.data.nwells;
            x=hDV.data.fault.xf(j); y=hDV.data.fault.yf(j);
            ppOnFault(j) = pfieldcalc(hDV,x,y,wells);
            fraction=(2/3)+(j/(3.*nfaults));
            waitbar(fraction,waitBarHandleWells,['Calculating ',num2str(100*fraction,'%3.0f'), '% ...']) ;
        end
        
        hDV.plotdata.pint.ppf=ppOnFault;
        dpth=hDV.data.stress.vals(5);
        
        sig = hDV.data.stress.vals(1:3)*dpth ; 
        
        pp0 = hDV.data.stress.vals(6)*dpth ;
        nu = hDV.data.reservoir.vals(5)  ;
        biot = 1 ; %biot Coefficient
        
        if hDV.data.stress.aphi.use
            inputCell = {sig',0.00,pp0,thf,dips,SHdir,ppOnFault,muf,biot,nu,hDV.data.stress.aphi.vals(1)} ;% add aphi
            % note that in this case, horizontal stress(es) in sig are
            % ignored, but entered - it is still a 3 element vector, but
            % the sHmax and/or shmin will be calculated using A-Phi or
            % Modified A-PHi model
        else
            inputCell = {sig',0.00,pp0,thf,dips,SHdir,ppOnFault,muf,biot,nu};
        end
        [~ , hDV.plotdata.resultsH.outs , hDV.plotdata.resultsH.C1 , hDV.plotdata.resultsH.C2 , hDV.plotdata.resultsH.C3 , ...
            hDV.plotdata.resultsH.sig_fault , hDV.plotdata.resultsH.tau_fault]=mohrs_3D(inputCell,hDV) ;
        
        if ishandle(waitBarHandleWells)
            close(waitBarHandleWells) ;
        end
        else % importing hydrology % 0 does hydrology calculation internally, 1 imports XYZPT model like from MODFLOW
           
            % interpolate pressure on faults
            ts = get(hDV.hdsldr(1),'value');
            numbersImportedHydrology=hDV.data.reservoir.numbersImportedHydrology;
            
            thisYearData=numbersImportedHydrology(numbersImportedHydrology(:,4)==ts,:);
            thisYearX=thisYearData(:,1);
            thisYearY=thisYearData(:,2);
            thisYearPSI=thisYearData(:,3);
            
            % interpolate imported hydrology onto fault
            ppOnFault=griddata(thisYearX,thisYearY,thisYearPSI,hDV.data.fault.xf,hDV.data.fault.yf) ;
            hDV.plotdata.pint.ppf= ppOnFault;
            
            %%% interpolate imported hydrology on to field surface
            xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
            ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
            npts=50 ;
            [hDV.plotdata.pflot.Xgrid, hDV.plotdata.pflot.Ygrid] = meshgrid(linspace(xmin,xmax,npts) , linspace(ymin,ymax,npts)) ;
            hDV.plotdata.pflot.Zgrid = griddata(thisYearX,thisYearY,thisYearPSI,hDV.plotdata.pflot.Xgrid,hDV.plotdata.pflot.Ygrid) ; % in PSI
            
            % below is all as above - could bring outside if statement
            hDV.plotdata.pint.Xgrid = hDV.plotdata.pflot.Xgrid;
            hDV.plotdata.pint.Ygrid = hDV.plotdata.pflot.Ygrid;
            hDV.plotdata.pint.Zgrid = hDV.plotdata.pflot.Zgrid;
            
        
            %fault failure calculation, Mohrs circle
            thf=hDV.data.fault.thf ;   %fault strikes
            dips=hDV.data.fault.dipf; %fault dips
            SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
            muf=hDV.data.fault.muf; %fault mu
            
            dpth=hDV.data.stress.vals(5);
            sig = hDV.data.stress.vals(1:3)*dpth ; 
            %
            pp0 = hDV.data.stress.vals(6)*dpth ;
            nu = hDV.data.reservoir.vals(5)  ;
            biot = 1 ; %biot Coefficient
            
            if hDV.data.stress.aphi.use
                inputCell = {sig',0.00,pp0,thf,dips,SHdir,ppOnFault,muf,biot,nu,hDV.data.stress.aphi.vals(1)} ;% add aphi
                % note that in this case, horizontal stress(es) in sig are
                % ignored, but entered - it is still a 3 element vector, but
                % the sHmax and/or shmin will be calculated using A-Phi or
                % Modified A-PHi model
            else
                inputCell = {sig',0.00,pp0,thf,dips,SHdir,ppOnFault,muf,biot,nu};
            end
        
            [~ , hDV.plotdata.resultsH.outs , hDV.plotdata.resultsH.C1 , hDV.plotdata.resultsH.C2 , hDV.plotdata.resultsH.C3 , ...
                hDV.plotdata.resultsH.sig_fault , hDV.plotdata.resultsH.tau_fault]=mohrs_3D(inputCell,hDV) ;
            % above is all as above - could bring outside if statement  
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'GEOMECHANICS' , 'PROB. GEOMECH', 'MODEL INPUTS'}
        
        %fault failure calculation, Mohrs circle
        thf=hDV.data.fault.thf ;   %fault strikes
        dips=hDV.data.fault.dipf; %fault dips
        SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
        
        muf=hDV.data.fault.muf; %fault friction coefficients
        dpth=hDV.data.stress.vals(5); %calculation depth
        sig = hDV.data.stress.vals(1:3)*dpth ; 
        
        
        setstressregtext(hDV);
        
        pp0 = hDV.data.stress.vals(6)*dpth ; %Native pore pressure
        nu = hDV.data.reservoir.vals(5)  ; %Poissons Ratio
        biot = 1 ; %biot Coefficient hard coded
        
        if hDV.data.stress.aphi.use
            inputCell = {sig',0.00,pp0,thf,dips,SHdir,0*thf,muf,biot,nu,hDV.data.stress.aphi.vals(1)} ;% add aphi
            % note that in this case, horizontal stress(es) in sig are
            % ignored, but entered - it is still a 3 element vector, but
            % the sHmax and/or shmin will be calculated using A-Phi or
            % Modified A-PHi model
        else
            inputCell ={sig',0.00,pp0,thf,dips,SHdir,0*thf,muf,biot,nu};
        end
            
        [~ , hDV.plotdata.results.outs , hDV.plotdata.results.C1 , hDV.plotdata.results.C2 , hDV.plotdata.results.C3 , ...
            hDV.plotdata.results.sig_fault , hDV.plotdata.results.tau_fault]=mohrs_3D(inputCell,hDV) ;
        
        %obtain results for composite plot
        Ngrid=50;
        strs=linspace(0,360,Ngrid) ; dips=linspace(0,90,Ngrid); [st_grid,dp_grid] = meshgrid(strs,dips);
        [xx,yy] = pol2cart((st_grid)*pi/180,dp_grid*pi/180/(pi/2)); %convert to radians and scale dip of pi/2 to 1
        ppval=0*xx;
        
        for n=1:1:Ngrid
           
            if hDV.data.stress.aphi.use
                % add aphi;  see above if block for inputCell for more details
                indatacell = {sig',0.00,pp0,st_grid(:,n),dp_grid(:,n),360-SHdir,0*ones(Ngrid,1),mean(muf)*ones(Ngrid,1),biot,nu,hDV.data.stress.aphi.vals(1)} ; %!!! Need to understand why you need 360-SHdir to make this work
            else
                indatacell = {sig',0.00,pp0,st_grid(:,n),dp_grid(:,n),360-SHdir,0*ones(Ngrid,1),mean(muf)*ones(Ngrid,1),biot,nu} ; %!!! Need to understand why you need 360-SHdir to make this work
            end
            [~,outs,~,~,~,~,~]=mohrs_3D(indatacell,hDV);
            ppval(:,n) = outs.ppfail;
        end
        
        ppval(ppval<0)=0; %no negatives in composite plot
        hDV.plotdata.results.compx=xx;
        hDV.plotdata.results.compy=yy;
        hDV.plotdata.results.compz=ppval;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end % switch 


end