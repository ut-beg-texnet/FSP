% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Refresh each plot panel

function refreshplotdata(hDV,name)

if hDV.plotdata.printFunctionName
    disp(['running refreshplotdata for ',name])
end

switch name
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'MODEL INPUTS'
        refreshFaultSelectorList(hDV,name) % refresh multi fault selector
        
        hp = hDV.plotdata.inputMap ;
        
        %setup axes for faults plot
        xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
        ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
        axis(hp.ax1,[xmin xmax ymin ymax]);
        axis(hp.ax1,'equal');
        
        %stress regime text
        setstressregtext(hDV);
        
        %faults data
        lenf=hDV.data.fault.lenf ;
        xf = hDV.data.fault.xf ; yf = hDV.data.fault.yf ;
        thf_map=hDV.data.fault.thf ;
        thf_cart=90-thf_map ; %convert to cartesian
        fstart = xf+1i*yf - lenf/2.*exp(1i*thf_cart*pi/180) ; fend   = xf+1i*yf + lenf/2.*exp(1i*thf_cart*pi/180) ;
        nfaults = hDV.data.fault.vals(1) ;
        
        %update the plot  for each fault
        for j=1:1:nfaults
            %fault plot
            xp = real([fstart(j) fend(j)]) ; yp = imag([fstart(j) fend(j)]) ;
            set(hDV.plotdata.flinesgeoModelInputs(j),'xdata',xp ,'ydata' , yp  ) ;
            set(hDV.plotdata.flinesgeobackModelInputs(j),'xdata',xp ,'ydata' , yp  ) ;
            set(hDV.plotdata.flinesgeoModelInputs(j), 'color',0.3.*[1,1,1]);
        end
        %dont plot remaining empty fault placeholders
        set(hDV.plotdata.flinesgeoModelInputs(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.flinesgeobackModelInputs(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        
        
        %refresh well labels annotation and X positions
        wellColors=lines(hDV.data.nwells);
        
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
            xsw=hDV.data.realWellData.XEasting;
            ysw=hDV.data.realWellData.YNorthing;
        else
            xsw=hDV.data.inject.vals(1:hDV.data.nwells,1);
            ysw=hDV.data.inject.vals(1:hDV.data.nwells,2);
        end
        set(hp.hwlocs,'xdata',xsw,'ydata',ysw);
        
        %anotations for wells
        lt=get(hp.ax1,'xlim') ;
        spc=0; %.005*(lt(2)-lt(1)) ; % space between point and text location, 
        % specifying spc as 0 and making space below before wstr means
        % label stays on zoomin
        for k=1:1:hDV.data.nwells
            pos = [spc+xsw(k) , spc+ysw(k)] ;
            
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                wstr = hDV.data.realWellData.wellNames{k};
            else %if hDV.data.nwells_max<=25
                %                   wstr = char(65+k-1);
                %  else % if more than 25 wells possible, number wells, don't letter them
                wstr = num2str(k);
            end
            
            set( hDV.plotdata.inputMap.hwellLabel(k),'Position',pos,'color',wellColors(k,:),'string',[' ',wstr]) ;
        end
        set( hDV.plotdata.inputMap.hwellLabel((1+hDV.data.nwells):1:hDV.data.nwells_max),'Position',[NaN NaN]) ;
        
        
        % wells through time lines
        startyear=get(hDV.hdsldr(1),'min');
        endyear=get(hDV.hdsldr(1),'max'); % this is hardwired in slider bar
        deltaPlot=1; % jump in years
        mxar = zeros(hDV.data.nwells,1);
        
        if hDV.data.reservoir.importHydrology==0
            for k=1:1:hDV.data.nwells
                if  isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use % variable or constant data?
                    datenumBarrelsPerDayThisWell=hDV.data.realWellData.datenumBarrelsPerDay{k,1};
                else
                    datenumBarrelsPerDayThisWell=hDV.data.inject.datenumBarrelsPerDay{k,1};
                end

                [tt,QQ]=stairs(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2));
                % plot(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2),'rx') % data points as red xs

                % do this adjustment because value represents rate up to that point, not
                % rate after that point as stair function assumes.
                % see plottingInjectionStairData.m in test code and SpreadsheetStrings2WellData
                QQ(1:2)=[];
                QQ(end+1:end+2,1)=QQ(end);
                % these points bring the line down to the X axis at the beginning and end
                ys=[0;QQ;0];
                tt=[tt(1);tt;tt(end)];
                % plot(tt,ys,'b') % injection pattern in blue
                xs=decyear(tt);

                mxar(k) = max(ys) ; %store max to set yaxis limits
                set(hDV.plotdata.inputMap.hwlocs2(k),'xdata',xs,'ydata',ys,'color',wellColors(k,:))
            end
            set(hDV.plotdata.inputMap.hwlocs2(hDV.data.nwells+1:1:hDV.data.nwells_max),'xdata',NaN,'ydata',NaN)  ;

            % axis2 limits
            set(hDV.plotdata.inputMap.ax2,'xlim',[ startyear,endyear]);
            ymaxval = max(mxar) ; if ymaxval==0, ymaxval=100; end
            set(hDV.plotdata.inputMap.ax2,'ylim',[ 0,1.2*ymaxval]);

            %vertical year bar set
            set(hDV.plotdata.inputMap.hWellTimeBar,'xdata',[1;1].*get(hDV.hdsldr(1),'value'),'ydata',get(hDV.plotdata.inputMap.ax2,'ylim'));


            %setup well number text
            xstr = cell(hDV.data.nwells+2,1) ; xstr{1}='All' ;  xstr{2}='NONE' ;
            for k=1:1:hDV.data.nwells
                if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                    xstr{k+2}=['Well #',hDV.data.realWellData.wellNames{k}];
                else
                    xstr{k+2}=['Well #',num2str(k)];
                end
            end
            set(hDV.plotdata.inputMap.popWells,'string',xstr) ;
        
        else
            % Hydrology model improted; no wells used
            %
            xstr = cell(1,1); xstr{1}='NONE';
            set(hDV.plotdata.inputMap.popWells,'string',xstr) ;
        end
        
        %stress direction arrow
        setuparrowstressdir(hDV);
        
        %        make single well plots visible?
        if hDV.data.reservoir.importHydrology==0 % if using well data, not imported hydrologic model
            set(hDV.plotdata.inputMap.ax2,'visible','on')
        elseif hDV.data.reservoir.importHydrology==1 % if imported hydrologic model, don't show single well solutions
            set(hDV.plotdata.inputMap.ax2,'visible','off')
            set(hDV.plotdata.inputMap.hwlocs2(1:hDV.data.nwells_max),'Xdata',NaN) ;
            set(hDV.plotdata.inputMap.hwlocs2(1:hDV.data.nwells_max),'Ydata',NaN) ;
            set(hp.hwlocs,'xdata',NaN,'ydata',NaN);
            set( hDV.plotdata.inputMap.hwellLabel(1:1:hDV.data.nwells_max),'Position',[NaN NaN])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'HYDROLOGY'
        
        %run the geomechanics calculations in case they weren't run
        calcengine(hDV,'GEOMECHANICS');
        
               
        hp = hDV.plotdata.pflot ; % hydrology contour plot
        colormap(hp.ax,'default') ;
        %%% refresh pressure field surface
        xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
        ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
        
        set(hp.hsf,'xdata',hp.Xgrid, ...
            'ydata',hp.Ygrid, ...
            'zdata',hp.Zgrid) ;
        
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
            xs=hDV.data.realWellData.XEasting(1:hDV.data.nwells);
            ys=hDV.data.realWellData.YNorthing(1:hDV.data.nwells);
        else
            xs=hDV.data.inject.vals(1:hDV.data.nwells,1);
            ys=hDV.data.inject.vals(1:hDV.data.nwells,2);
        end
        
        zs=0*xs+max(hp.Zgrid(:));
        set(hp.hwlocs,'xdata',xs,'ydata',ys,'zdata',zs) ;
        
        %refresh well labels annotation
        lt=get(hp.ax,'xlim') ;
        spc=.01*(lt(2)-lt(1)) ;
        for k=1:1:hDV.data.nwells
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                xp=hDV.data.realWellData.XEasting(k);
                yp=hDV.data.realWellData.YNorthing(k);
            else
                xp=hDV.data.inject.vals(k,1) ;  yp=hDV.data.inject.vals(k,2)  ;
            end
            pos = [spc+xp , spc+yp  ,  zs(1)] ;
            
            
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                wstr = hDV.data.realWellData.wellNames{k};
            else
                wstr = num2str(k);  % char(65+k-1);
            end
            set(hp.htx(k),'Position',pos,'string',wstr) ;
        end
        set(hp.htx(1+hDV.data.nwells:1:hDV.data.nwells_max),'Position',[NaN NaN NaN]) ;
        axis(hp.ax,'equal','tight');
        axis(hp.ax,[xmin xmax ymin ymax]);
        
        %%% refresh individual well curves
        for k=1:1:hDV.data.nwells
            set(hp.pl(k),'Xdata',hp.x(k,:),'Ydata',hp.y(k,:));
        end
        set(hp.pl(hDV.data.nwells+1:hDV.data.nwells_max),'Xdata',NaN) ;
        set(hp.pl(hDV.data.nwells+1:hDV.data.nwells_max),'Ydata',NaN) ;
        
        
        %%% Mohrs circle moving
        dpth = hDV.data.stress.vals(5); sig=hDV.data.stress.vals(1:3)*dpth ;
        nfaults = hDV.data.fault.vals(1) ;
        xlim(hp.ax3,[0 1.2*max(sig)]); ylim(hp.ax3,[0 (max(sig)-min(sig))]) ;
        
        cv=hDV.plotdata.results.outs.ppfail;
        %mincfc=0 ; maxcfc=min(dpth*hDV.data.stress.vals(1:3)) ;
        
            for j=1:1:nfaults
                if cv(j)<=0
                    clr = [1 0 0] ;
                else
                    clr = .5*[1 1 1] ;
                end
                
                set(hp.mc1(j),'xdata',real(hDV.plotdata.resultsH.C1(j,:)) ,'ydata',imag(hDV.plotdata.resultsH.C1(j,:)),'color',clr,'visible','off');
                set(hp.mc2(j),'xdata',real(hDV.plotdata.resultsH.C2(j,:)) ,'ydata',imag(hDV.plotdata.resultsH.C2(j,:)),'color',clr,'visible','off');
                set(hp.mc3(j),'xdata',real(hDV.plotdata.resultsH.C3(j,:)) ,'ydata',imag(hDV.plotdata.resultsH.C3(j,:)),'color',clr,'visible','off');
            end
            
            for j=1:1:nfaults
                cl = getcolor(flipud(hDV.cmapGYR),cv(j),hDV.plotdata.mincfc,hDV.plotdata.maxcfc) ;
                set(hp.mflt(j),'xdata',hDV.plotdata.resultsH.sig_fault(j), 'ydata' , hDV.plotdata.resultsH.tau_fault(j));
                set(hp.mflt(j),'markerfacecolor',cl,'visible','off') ;
            end
 
        set(hp.mc1(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;  %dont plot remaining empty fault placeholders
        set(hp.mc2(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hp.mc3(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hp.mflt(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        uistack(hp.mflt(:),'top') ;
        x=linspace(0,max(sig)) ; mu =  hDV.data.fault.vals(2) ;
        set(hp.mfln , 'xdata' , x ,  'ydata' , x*mu) ;
                   
        refreshFaultSelectorList(hDV,'HYDROLOGY') % refresh multi fault selector, this makes visible selected faults on mohr circle
        callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV,'HYDROLOGY')
        
        %setup well number text
        xstr = cell(hDV.data.nwells+2,1) ; xstr{1}='All' ;  xstr{2}='None' ;
        for k=1:1:hDV.data.nwells
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                xstr{k+2}=['Well #',hDV.data.realWellData.wellNames{k}];
            else
                xstr{k+2}=['Well #',num2str(k)];
            end
        end
        set(hDV.plotdata.pflot.popWells,'string',xstr) ;
        
        % make single well plots visible?
        if hDV.data.reservoir.importHydrology==0 % if using well data, not imported hydrologic model
            set(hDV.plotdata.pflot.ax2,'visible','on')
        elseif hDV.data.reservoir.importHydrology==1 % if imported hydrologic model, don't show single well solutions
            set(hDV.plotdata.pflot.ax2,'visible','off')
            set(hp.pl(1:hDV.data.nwells_max),'Xdata',NaN) ;
            set(hp.pl(1:hDV.data.nwells_max),'Ydata',NaN) ;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'GEOMECHANICS'
        hp = hDV.plotdata.pffot ; %select panel
        colormap(hDV.plotdata.pffot.ax3,flipud(hDV.cmapGYR)) ;
        %setup axes for faults map plot
        xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
        ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
        axis(hp.ax,[xmin xmax ymin ymax]);
        
        %faults data
        lenf=hDV.data.fault.lenf ;
        xf = hDV.data.fault.xf ;
        yf = hDV.data.fault.yf ;
        thf_map=hDV.data.fault.thf ; dipf=hDV.data.fault.dipf ;
        thf_cart=90-thf_map ; %convert to cartesian
        fstart = xf+1i*yf - lenf/2.*exp(1i*thf_cart*pi/180) ; fend   = xf+1i*yf + lenf/2.*exp(1i*thf_cart*pi/180) ;
        nfaults = hDV.data.fault.vals(1) ;
        

        %color range and value to use
        dpth = hDV.data.stress.vals(5);
        cv=hDV.plotdata.results.outs.ppfail;
        hDV.plotdata.mincfc = str2double(get(hDV.plotdata.pffot.cmintxt,'string'));
        hDV.plotdata.maxcfc = str2double(get(hDV.plotdata.pffot.cmaxtxt,'string'));
        
        if isnan(hDV.plotdata.mincfc) || isnan(hDV.plotdata.maxcfc)
            hDV.plotdata.mincfc=0 ; hDV.plotdata.maxcfc=min(hDV.data.stress.vals(1:3)*dpth);
            set(hDV.plotdata.pffot.cmintxt,'string',num2str(hDV.plotdata.mincfc,'%5.2f'));
            set(hDV.plotdata.pffot.cmaxtxt,'string',num2str(hDV.plotdata.maxcfc));
        end
        
        %update the plot and stereonet data for each fault
        for j=1:1:nfaults
            %fault plot
            xp = real([fstart(j) fend(j)]) ; yp = imag([fstart(j) fend(j)]) ;
            set(hDV.plotdata.flinesgeo(j),'xdata',xp ,'ydata' , yp  ) ;
            set(hDV.plotdata.flinesgeoback(j),'xdata',xp ,'ydata' , yp  ) ;
            cl = getcolor(flipud(hDV.cmapGYR),cv(j),hDV.plotdata.mincfc,hDV.plotdata.maxcfc) ;
            set(hDV.plotdata.flinesgeo(j), 'color', cl);
            %stereonet
            [Theta, Rho, pTh, pRh]=stereonet(dipf(j),thf_cart(j));
            [X,Y] = pol2cart(Theta,Rho);
            [Xp,Yp] = pol2cart(pTh,pRh);% poles
            set(hDV.plotdata.snet(j),'xdata',X ,'ydata' , Y ) ;
            set(hDV.plotdata.snetpoles(j),'xdata',Xp ,'ydata' , Yp ) ;
            set(hDV.plotdata.snet(j),'color',cl) ;
            set(hDV.plotdata.snetpoles(j),'markerfacecolor',cl,'markeredgecolor',cl) ;
            %mohrs circle
            set(hDV.plotdata.pffot.mflt(j),'xdata',hDV.plotdata.results.sig_fault(j), 'ydata' , hDV.plotdata.results.tau_fault(j));
            set(hDV.plotdata.pffot.mflt(j),'markerfacecolor',cl) ;
            
        end
        
        %equal axes
        axis(hp.ax,'equal');
        
        %dont plot remaining empty fault placeholders
        set(hDV.plotdata.flinesgeo(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.flinesgeoback(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.snetpoles(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.snet(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.pffot.mflt(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        uistack(hDV.plotdata.pffot.mflt(:),'top') ;
        
        %max hor stress direction display
        %         th=90-hDV.data.stress.vals(4);
        %         [X,Y] = pol2cart(th*pi/180*[1 1],[-1 1]);
        [X,Y] = setuparrowstressdirStereonet(hDV);
        set(hDV.plotdata.sHmaxdir,'xdata',X ,'ydata' , Y ) %max hor direction
        set(hDV.plotdata.sHmaxdir,'color',[.5 .5 .5],'linewidth',2,'linestyle','-') ;
        
        %composite poles plot
        set(hDV.plotdata.snetcomp,'xdata',hDV.plotdata.results.compx,...
            'ydata',hDV.plotdata.results.compy,...
            'zdata',hDV.plotdata.results.compz);
        caxis(hDV.plotdata.pffot.ax3,[hDV.plotdata.mincfc,hDV.plotdata.maxcfc]) ;
        
        %do not show text on plots on panel change (slows things down
        %otherwise)
        set(hDV.plotdata.flinesgeotxt(:),'visible', 'off') ;
        set(hDV.plotdata.pffot.pop(1),'value',1) ;
        
        %stereonet labels
        set(findall(gcf, 'String', '  1','-or','String', '  0.8'),'String', '   '); %remove the radial labels
        set(findall(gcf, 'String', '  0.6','-or','String', '  0.4'),'String', '   '); %remove the radial labels
        set(findall(gcf, 'String', '  0.2','-or','String', '  0.5'),'String', '   '); %remove the radial labels
        
        
        set(findall(gcf, 'String', '0'),'String', ' E','clipping','on');
        set(findall(gcf, 'String', '30'),'String', ' 60','clipping','on');
        set(findall(gcf, 'String', '60'),'String', ' 30','clipping','on');
        set(findall(gcf, 'String', '90'),'String', ' N','clipping','on');
        set(findall(gcf, 'String', '120'),'String', ' 330','clipping','on');
        set(findall(gcf, 'String', '150'),'String', ' 300','clipping','on');
        set(findall(gcf, 'String', '180'),'String', ' W','clipping','on');
        set(findall(gcf, 'String', '210'),'String', ' 240','clipping','on');
        set(findall(gcf, 'String', '240'),'String', ' 210','clipping','on');
        set(findall(gcf, 'String', '270'),'String', ' S','clipping','on');
        set(findall(gcf, 'String', '300'),'String', ' 150','clipping','on');
        set(findall(gcf, 'String', '330'),'String', ' 120','clipping','on');
        
        %mohrs circle
        sig=hDV.data.stress.vals(1:3)*dpth ; pp=hDV.data.stress.vals(6)*dpth;
        xlim(hp.ax2,[0 1.2*(max(sig)-pp)]); ylim(hp.ax2,[0 (max(sig)-min(sig))]) ;
        
        clr = .5*[1 1 1] ; kc=1 ;
        set(hDV.plotdata.pffot.mc1,'xdata',real(hDV.plotdata.results.C1(kc,:)) ,'ydata',imag(hDV.plotdata.results.C1(kc,:)),'color',clr );
        set(hDV.plotdata.pffot.mc2,'xdata',real(hDV.plotdata.results.C2(kc,:)) ,'ydata',imag(hDV.plotdata.results.C2(kc,:)),'color',clr );
        set(hDV.plotdata.pffot.mc3,'xdata',real(hDV.plotdata.results.C3(kc,:)) ,'ydata',imag(hDV.plotdata.results.C3(kc,:)),'color',clr  );
        x=linspace(0,max(sig)) ; mu =  hDV.data.fault.vals(2) ;
        set(hDV.plotdata.pffot.mfln , 'xdata' , x ,  'ydata' , x*mu) ;
        
        %            sigma labels on mohr circle
        sigmaLabelYs=[0] ;%[1;1;1]*yMax*0.9
        sigmaLabelXs=[sig-pp];
        set(hDV.plotdata.pffot.SigmaHhVLabel1,'position',[sigmaLabelXs(1) , sigmaLabelYs])
        set(hDV.plotdata.pffot.SigmaHhVLabel2,'position',[sigmaLabelXs(2) , sigmaLabelYs])
        set(hDV.plotdata.pffot.SigmaHhVLabel3,'position',[sigmaLabelXs(3) , sigmaLabelYs])
        %
        
        %refresh the color bar
        callcbartxt([],[],hDV);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'PROB. GEOMECH'
        %check if the button for prob calc is red
        if all(get(hDV.plotdata.pprob.hbutMC,'backgroundcolor')==hDV.colors.red)
            %run the prob calc
            callrunmc(hDV.plotdata.pprob.hbutMC ,[],hDV);
        end
        refreshFaultSelectorList(hDV,'PROB. GEOMECH') % multi fault selector
        callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV)
        updateinputbar(hDV); %update the input variability bar plot
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case  'PROB. HYDRO'
        hp = hDV.plotdata.hydprob;
        nfaults = hDV.data.fault.vals(1) ;
        
        % run the PROB. GEOMECH if it's red
        if all(get(hDV.plotdata.pprob.hbutMC,'backgroundcolor')==hDV.colors.red)
            callrunmc(hDV.plotdata.pprob.hbutMC ,[],hDV);
        end
        
        if ~isfield(hDV.data.probHydrology,'distirbutionTxt') % if this field doesnt exist, calldatHydroMc to make you enter parameters
            calldatHydroMc([],[],hDV)
            return
        end
        
        % run probabilistic hydrology if button is red, or if plotted year
        % isn't slider year
        if all(get(hDV.plotdata.hydprob.hbutMC,'backgroundcolor')==hDV.colors.red)  ||  hDV.plotdata.hydprob.yearDataShown~=get(hDV.hdsldr(3),'val')
            callRunHydro_mc([],[],hDV)
        end
        
        [hDV.plotdata.pint.fsp]=calcFaultFSP(hDV); % calculate FSP for all faults
        
        % populate dropdown menu with colred fault background
        if ~ isfield(hDV.plotdata,'minint')
            hDV.plotdata.minint=0;
            hDV.plotdata.maxint=1;
        end
        
        refreshFaultSelectorList(hDV,'PROB. HYDRO') % multi fault selector
        callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV)
        
        % refresh background dashed CDF lines
        for jj=1:nfaults
            cl = getcolor(hDV.cmapGYR,hDV.plotdata.pint.fsp(jj),hDV.plotdata.minint,hDV.plotdata.maxint) ; % find RGB color
            set(hDV.plotdata.greyCDFBackground(jj),'xdata',get(hDV.plotdata.flinesprob(jj),'xdata'),'ydata',get(hDV.plotdata.flinesprob(jj),'ydata'),'color',cl)
        end
        set(hDV.plotdata.greyCDFBackground(nfaults+1:hDV.data.NFAULTSMAX),'xdata',NaN,'ydata',NaN)
        
        %         set(hDV.plotdata.hydprob.ax2,'xlim',
        
        
        
    case 'INTEGRATED'
        hp = hDV.plotdata.pint ; %select panel
        colormap(hp.ax,'default') ;
        %setup axes for faults map plot
        xmin=hDV.data.adv.vals(1); xmax=hDV.data.adv.vals(2);
        ymin=hDV.data.adv.vals(3); ymax=hDV.data.adv.vals(4);
        
        %refresh surface
        set(hp.hsf,'xdata',hp.Xgrid, ...
            'ydata',hp.Ygrid, ...
            'zdata',hp.Zgrid) ;
        
        %run the geomechanics calculations in case they weren't run
        calcengine(hDV,'GEOMECHANICS');
        
        %run the PROB. GEOMECH if its red
        if all(get(hDV.plotdata.pprob.hbutMC,'backgroundcolor')==hDV.colors.red)
            callrunmc(hDV.plotdata.pprob.hbutMC ,[],hDV);
        end
        
        % run probabilistic hydrology if button is red, or if plotted year
        % isn't slider year
        if all(get(hDV.plotdata.hydprob.hbutMC,'backgroundcolor')==hDV.colors.red)  || hDV.plotdata.hydprob.yearDataShown~=get(hDV.hdsldr(2),'val')
            callRunHydro_mc([],[],hDV)
        end
        
        [hDV.plotdata.pint.fsp]=calcFaultFSP(hDV); % calculate FSP for all faults
        
        %faults data
        lenf=hDV.data.fault.lenf ;
        xf = hDV.data.fault.xf ;
        yf = hDV.data.fault.yf ;
        thf_map=hDV.data.fault.thf ; dipf=hDV.data.fault.dipf ;
        thf_cart=90-thf_map ; %convert to cartesian
        fstart = xf+1i*yf - lenf/2.*exp(1i*thf_cart*pi/180) ; fend   = xf+1i*yf + lenf/2.*exp(1i*thf_cart*pi/180) ;
        nfaults = hDV.data.fault.vals(1) ;
        
        %update the main map plot
        for j=1:1:nfaults
            %fault plot
            xp = real([fstart(j) fend(j)]) ; yp = imag([fstart(j) fend(j)]) ;
            set(hDV.plotdata.flinesintback(j),'xdata',xp ,'ydata' , yp  ) ;
            set(hDV.plotdata.flinesint(j),'xdata',xp ,'ydata' , yp  ) ;
        end
        
        % dont plot remaining empty fault placeholders
        set(hDV.plotdata.flinesintback(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        set(hDV.plotdata.flinesint(nfaults+1:end),'xdata', NaN ,'ydata' , NaN ) ;
        
        %equal axes
        axis(hp.ax,'equal','tight');
        axis(hp.ax,[xmin xmax ymin ymax]);
        
        %compute the pressure margin and fsp
        limits=hDV.plotdata.results.outs.ppfail(1:nfaults); %obtain limits
        if all(isnan(limits))
            hDV.plotdata.pint.ppm=NaN;
            errordlg('Data not current! Run PROB. GEOMECH tab to compute pressure margins');
        else
            ppm = limits-hDV.plotdata.pint.ppf ; ppm(ppm<0)=0;
            hDV.plotdata.pint.ppm=ppm;
            %             hDV.plotdata.pint.fsp=0*ppm;
            
            %             [hDV.plotdata.pint.fsp]=calcFaultFSP(hDV);
            
            %             for j=1:1:nfaults
            %                 pres = get(hDV.plotdata.flinesprob(j),'xdata');
            %                 prob = get(hDV.plotdata.flinesprob(j),'ydata');
            %                 [~,ix]=unique(pres); pres = pres(ix) ; prob = prob(ix) ;
            %
            %                 if length(pres)==1 %this is the deterministic case
            %                     if hDV.plotdata.pint.ppf(j)<pres
            %                         hDV.plotdata.pint.fsp(j)=0;
            %                     else
            %                         hDV.plotdata.pint.fsp(j)=1;
            %                     end
            %                 else
            %                     hDV.plotdata.pint.fsp(j) = interp1(pres,prob,hDV.plotdata.pint.ppf(j),'nearest','extrap');
            %                 end
            %             end
            
            
            %             set(hp.pop(2),'value',hDV.plotdata.curfault(1)+1);
            %run popup response again
            %             callintpop(hDV.plotdata.pint.pop(1),[],hDV);
        end
        
        %for color ranges
        hDV.plotdata.minint = str2double(get(hDV.plotdata.pint.cmintxt,'string'));
        hDV.plotdata.maxint = str2double(get(hDV.plotdata.pint.cmaxtxt,'string'));
        
        if isnan(hDV.plotdata.minint) || isnan(hDV.plotdata.maxint)
            hDV.plotdata.minint=0 ; hDV.plotdata.maxint=1;
            set(hDV.plotdata.pint.cmintxt,'string',num2str(hDV.plotdata.minint,'%5.2f'));
            set(hDV.plotdata.pint.cmaxtxt,'string',num2str(hDV.plotdata.maxint));
        end
        
        
        if hDV.data.reservoir.importHydrology==0
            if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use
                xsw=hDV.data.realWellData.XEasting;
                ysw=hDV.data.realWellData.YNorthing;
            else
                xsw=hDV.data.inject.vals(1:hDV.data.nwells,1);
                ysw=hDV.data.inject.vals(1:hDV.data.nwells,2);
            end
            set(hDV.plotdata.pint.hwlocs,'xdata',xsw,'ydata',ysw);
        else % pressure model imported, don't show wells
            set(hDV.plotdata.pint.hwlocs,'xdata',NaN,'ydata',NaN);
        end
        
        
        %update colors of faults
        callintcbartxt([],[],hDV); % this also updates faults selected
        
        %year axes on right side plots
        tlims= [get(hDV.hdsldr(2),'min')-0.001 , get(hDV.hdsldr(2),'max')+0.001];
        xlim(hDV.plotdata.pint.ax2,tlims);
        xlim(hDV.plotdata.pint.ax4,tlims)
        %green time bars refresh
        set(hDV.plotdata.pint.timebar1,'xdata',[1;1].*get(hDV.hdsldr(1),'value'),'ydata',get(hDV.plotdata.pint.ax2,'ylim'));
        set(hDV.plotdata.pint.timebar2,'xdata',[1;1].*get(hDV.hdsldr(1),'value'),'ydata',get(hDV.plotdata.pint.ax4,'ylim'));
        
end


end