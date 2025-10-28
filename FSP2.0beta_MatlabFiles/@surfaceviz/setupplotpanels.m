% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Setup the different GUI plot panels. See surfaceviz.m
% modified by Rall


function setupplotpanels(hDV)

for k=1:1:length(hDV.tabnames)
    pl=hDV.pMain(k); %which panel am I in
    switch hDV.tabnames{k}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'MODEL INPUTS'
            
            %show stress regime in text
            hDV.plotdata.inputMap.srtxt=uicontrol('parent',pl,'style','text','string','Stress Regime: Strike Slip Faulting'...
                ,'position',[.075 .85 .35 .05],'fontsize',14,'fontunits','normalized','backgroundcolor',[1 1 1]);
            
            %main plot axes to show faults
            hDV.plotdata.inputMap.ax1=axes('parent',pl,'position',[.07  .3  .38  .45],'fontsize',hDV.ftsz) ;
            hold(hDV.plotdata.inputMap.ax1,'on') ;
            xlabel(hDV.plotdata.inputMap.ax1,'x easting [km]','fontsize',hDV.ftsz+2) ; % why don't these stay?
            ylabel(hDV.plotdata.inputMap.ax1,'y northing [km]','fontsize',hDV.ftsz+2) ;
            grid(hDV.plotdata.inputMap.ax1,'on');
            
            
            %plotting black faults placeholder
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.flinesgeobackModelInputs(j)=plot(0,0,'linewidth',2,'color',hDV.colors.darkgrey,'parent',hDV.plotdata.inputMap.ax1) ;
                hDV.plotdata.flinesgeoModelInputs(j)=plot(0,0,'linewidth',3,'color','k','parent',hDV.plotdata.inputMap.ax1,'visible','off') ;
                hDV.plotdata.flinesgeotxtModelInputs(j)=text(0,0,' ');
                
            end
            set(hDV.plotdata.flinesgeotxtModelInputs(:),'visible','off') ;
            
            %text annotation of well locations
            for j=1:1:hDV.data.nwells_max
                hDV.plotdata.inputMap.hwellLabel(j)=text(NaN,NaN,0,char(65+j-1),'parent',hDV.plotdata.inputMap.ax1,'color',[1 1 1],'fontsize',20,'FontWeight','bold') ;
            end

            % X for each well
            hDV.plotdata.inputMap.hwlocs=plot(-1,-1,'ks','linewidth',2,'markersize',8,'parent',hDV.plotdata.inputMap.ax1) ;
            
            % clipping so they stay in axis when zooming
            set(hDV.plotdata.inputMap.hwellLabel(:),'clipping','on')     

            
            % SH azimuth plot
            % still need to fix aspect ratio
            hDV.plotdata.inputMap.ax3=axes('parent',pl,'position',[.07  .15  .1  .1],'fontsize',hDV.ftsz) ;
            set(hDV.plotdata.inputMap.ax3,'xtick',[],'ytick',[],'box','off','xlim',[-1,1],'ylim',[-1,1]) ;
            hDV.plotdata.inputMap.hPlotArrowSH=plot(NaN,NaN,'-','linewidth',2,'color','k');
            hold(hDV.plotdata.inputMap.ax3,'on');
            
            %dropdown menu control for selecting which wells to plot
            uicontrol('parent',pl,'style','text','string','Select Well:','position',[.55 .85 .2 .05],...
                'fontsize',14,'fontunits','normalized','backgroundcolor',[1 1 1]);
            hDV.plotdata.inputMap.popWells = uicontrol('parent',pl,'style','popup','string','All|Well #1|Well #2|Well #3','position',[.8 .85 .15 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callwellhighlight,hDV},'value',1,'backgroundcolor',[1 1 1]);
            
            %   wells through time plot
            hDV.plotdata.inputMap.ax2=axes('parent',pl,'position',[.57  .3  .38  .45],'fontsize',hDV.ftsz) ;
            %              set(hDV.plotdata.inputMap.ax2,'datetick','x')
            
            
            % wells through time line holders
            maxPlots = hDV.data.nwells_max ;
            hDV.plotdata.inputMap.hWellTimeBar=plot( hDV.plotdata.inputMap.ax2,NaN,NaN,'r--','linewidth',2,'visible','off') ; % ,'markersize',10
            hold(hDV.plotdata.inputMap.ax2,'on');
            for j=1:maxPlots
                hDV.plotdata.inputMap.hwlocs2(j)=plot(-1,-1,'-','color','k','linewidth',2) ; % ,'markersize',10
            end
            xlabel(hDV.plotdata.inputMap.ax2,'Time [years]','fontsize',hDV.ftsz+2) ;
            ylabel(hDV.plotdata.inputMap.ax2,'well rate [bbls/day]','fontsize',hDV.ftsz+2) ;
            grid(hDV.plotdata.inputMap.ax2,'on'); 
            
            %slider
            %GW=0.005; HE = 0.06 ; WC = 0.18 ; WM = 1 - (WC + 3 * GW); PI = [.005 0.005];
            %[hDV.hdsldr(1) hDV.hdsldr_txt(1)] = dateSldr(pl,PI,WM,HE,GW,hDV); %Slider #1 belongs to MODEL INPUTS panel
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'HYDROLOGY'
            
            mohrCircleLinewidth=1;
%             mohrCircleDotSize=10;

            %%% Surface plot of pressure field
            Xmax=10; Ymax=Xmax; npts=100 ; %grid spacing for plotting
            [hDV.plotdata.pflot.Xgrid , hDV.plotdata.pflot.Ygrid] = meshgrid(linspace(0,Xmax,npts) , linspace(0,Ymax,npts)) ;
            hDV.plotdata.pflot.Zgrid = 0*hDV.plotdata.pflot.Xgrid ;
            hDV.plotdata.pflot.ax=axes('parent',pl,'position',[.07  .25  .38  .5],'fontsize',hDV.ftsz) ;
            hDV.plotdata.pflot.hsf = surf(hDV.plotdata.pflot.Xgrid , hDV.plotdata.pflot.Ygrid , hDV.plotdata.pflot.Zgrid ,'parent',hDV.plotdata.pflot.ax);
            hold(hDV.plotdata.pflot.ax,'on') ; colormap(hDV.plotdata.pflot.ax,'default') ;
            shading(hDV.plotdata.pflot.ax,'interp');
            axis(hDV.plotdata.pflot.ax,'equal');
            axis(hDV.plotdata.pflot.ax,'tight');
            view(hDV.plotdata.pflot.ax,[0 90]) ;
            xlabel(hDV.plotdata.pflot.ax,'x easting [km]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pflot.ax,'y northing [km]','fontsize',hDV.ftsz) ;
            zlabel('Pressure front [psi]') ;
            h=colorbar ;
            ylabel(h,'Pressure front [psi]');
            


            hDV.plotdata.pflot.hwlocs=plot3(-1,-1,0,'kx','linewidth',2,'markersize',10) ;
            %text annotation of well locations
            for j=1:1:hDV.data.nwells_max
                hDV.plotdata.pflot.htx(j)=text(-1,-1,0,char(65+j-1),'parent',hDV.plotdata.pflot.ax,'color',[1 1 1],'fontsize',16,'FontWeight','normal') ;
            end
            cmax=1000 ; set(hDV.plotdata.pflot.ax,'clim',[0 cmax]) ; 
            
            %individual wells pressure curves
            hDV.plotdata.pflot.ax2=axes('parent',pl,'position',[.55  .48  .4  .42]) ;
            hold(hDV.plotdata.pflot.ax2,'all') ;
            maxPlots = hDV.data.nwells_max ;
            hDV.plotdata.pflot.x = zeros(maxPlots,2) ;  hDV.plotdata.pflot.y = zeros(maxPlots,2) ;
            hDV.plotdata.pflot.pl = zeros(maxPlots,1) ;
            
            for j=1:1:maxPlots
                hDV.plotdata.pflot.pl(j) = plot(hDV.plotdata.pflot.ax2,0,0,'linewidth',2) ;
            end
            
            xlabel('Distance [km]','fontsize',hDV.ftsz) ; ylabel('Pressure [psi]','fontsize',hDV.ftsz) ; grid(hDV.plotdata.pflot.ax2,'on');
             title('Single Well Radial Solutions','fontsize',hDV.ftsz,'fontunits','normalized') ;
            ymax=1000 ; ylim(hDV.plotdata.pflot.ax2,[0 ymax]) ;
            set(hDV.plotdata.pflot.ax2,'FontSize',hDV.ftsz);
            
            uicontrol('parent',pl,'style','text','string','Max plot DP [psi]:','position',[.05 .89 .2 .05],...
                'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1]);
            uicontrol('parent',pl,'style','edit','string',num2str(cmax),'position',[.25 .9 .1 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callmaxcpf,hDV});
            
            %dropdown menu control for selecting which wells to plot
            uicontrol('parent',pl,'style','text','string','Select Well: ','position',[0 .8 .22 .05],...
                'fontsize',14,'fontunits','normalized','backgroundcolor',[1 1 1]);
            hDV.plotdata.pflot.popWells = uicontrol('parent',pl,'style','popup','string','All|Well #1|Well #2|Well #3','position',[.22 .8 .15 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callwellhighlight2,hDV},'value',1,'backgroundcolor',[1 1 1],'horizontalalignment','right');
                        
            %MOHRs Circles moving
            hDV.plotdata.pflot.ax3=axes('parent',pl,'position',[.55 .05 .4 .3],'XColor','k','YColor','k','fontsize',hDV.ftsz,'color',hDV.colors.axsBgrnd) ;
            hold(hDV.plotdata.pflot.ax3,'on');
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.pflot.mc1(j)=plot(hDV.plotdata.pflot.ax3,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
                hDV.plotdata.pflot.mc2(j)=plot(hDV.plotdata.pflot.ax3,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
                hDV.plotdata.pflot.mc3(j)=plot(hDV.plotdata.pflot.ax3,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
            end
            

            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.pflot.mflt(j)=plot(hDV.plotdata.pflot.ax3,0,0,'ko','markerfacecolor','k','markersize',hDV.plotdata.mohrCircleDotSize);
            end

%             hDV.plotdata.pflot.chk = uicontrol('parent',pl,'style','checkbox','position',[.55 .34 .05 .05],...
%                 'backgroundcolor',[1 1 1],'callback',{@callchkmohr,hDV},'value',1);
%             uicontrol('parent',pl,'style','text','string','Show Mohr Circles','position',[.57 .35 .2 .03],...
%                 'fontsize',10,'fontunits','normalized','backgroundcolor',[1 1 1]);
%             
            hDV.plotdata.pflot.mfln =plot(0,0,'r','linewidth',2);
            axis(hDV.plotdata.pflot.ax3,'equal'); %view(hDV.plotdata.pmohr.ax,[0 90]);
            xlabel(hDV.plotdata.pflot.ax3,'\sigma effective [psi]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pflot.ax3,'\tau [psi]','fontsize',hDV.ftsz) ;
            grid(hDV.plotdata.pflot.ax3,'on');
            title(hDV.plotdata.pflot.ax3,'Mohrs Circles for All Faults','fontsize',12);
            
             %button to export CSV of CDF curves
             uicontrol('parent',pl,'style','push','string','Export Hydrology','position',[.05 .1 .16 .05],...
                'fontsize',hDV.ftsz,'fontunits','normalized','callback',{@callExportFSPHydrology_CSV,hDV},...
                'TooltipString','If you want to export a .csv file of X,Y,Pressure calculated in this year');
            
            
            %slider
            GW=0.005; HE = 0.06 ; WC = 0.18 ; WM = 1 - (WC + 3 * GW); PI = [.005 0.005];
            [hDV.hdsldr(1) hDV.hdsldr_txt(1)] = dateSldr(pl,PI,WM,HE,GW,hDV); %Slider #1 belongs to HYDROLOGY panel
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'GEOMECHANICS'
            mohrCircleLinewidth=2;
%             mohrCircleDotSize=10;
            %popup to choose what to show as labels no plot (CFF, SCU, PPF,
            %dip, etc.)
            hDV.plotdata.pffot.pop(1) = uicontrol('parent',pl,'style','popup','string',...
                'Choose Plot Labels|a) Fault Number|b) Fault Dip|c) Fault Strike|d) PPF: Pore Pressure to Slip|e) CFF: Coulomb Failure Function|f) SCU: Shear Capacity Utilization',...
                'position',[.05 .85 .3 .1],'fontsize',14,'fontunits','normalized','callback',{@callgeopop,hDV},'value',1);
            
            %help button
            hDV.plotdata.pffot.buthelp = uicontrol('parent',pl,'style','push','string','Help',...
                'position',[.37 .9 .08 .06],'fontsize',14,'fontunits','normalized','callback',{@callbuthelp,hDV,'help.jpg'},...
            'TooltipString','a schematic of what is being asessed in a Mohr Diagram');
  
            hDV.plotdata.curfault = 0 ; %current fault that is chosen for plotting, 0 is all
             
            %main plot axes to show faults
            hDV.plotdata.pffot.ax=axes('parent',pl,'position',[.07  .3  .38  .45],'fontsize',hDV.ftsz,'color',hDV.colors.axsBgrnd) ;
            
            hold(hDV.plotdata.pffot.ax,'on') ; 
            %plotting colored faults placeholder
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.flinesgeoback(j)=plot(0,0,'linewidth',2,'color',hDV.colors.darkgrey,'parent',hDV.plotdata.pffot.ax) ;
                hDV.plotdata.flinesgeo(j)=plot(0,0,'linewidth',3,'color','k','parent',hDV.plotdata.pffot.ax) ;
                hDV.plotdata.flinesgeotxt(j)=text(0,0,1,' ','parent',hDV.plotdata.pffot.ax,'fontsize',12,'clipping','on');
            end
            set(hDV.plotdata.flinesgeotxt(:),'visible','off') ;
            xlabel(hDV.plotdata.pffot.ax,'x easting [km]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pffot.ax,'y northing [km]','fontsize',hDV.ftsz) ;
            grid(hDV.plotdata.pffot.ax,'on');
            
            
            %controls for color bar range geomechanics
            hDV.plotdata.pffot.cmaxtxt = uicontrol('parent',pl,'style','edit','string','--','position',[.425 .07 .05 .03],...
                'fontsize',12,'fontunits','normalized','HorizontalAlignment','center','callback',{@callcbartxt,hDV});
            hDV.plotdata.pffot.cmintxt = uicontrol('parent',pl,'style','edit','string','--','position',[.025 .07 .05 .03],...
                'fontsize',12,'fontunits','normalized','HorizontalAlignment','center','callback',{@callcbartxt,hDV});
         
            % plot colorbar axis for geomech
            hDV.plotdata.pffot.axCBAR=axes('parent',pl,'position',[.05  .1  .4  .05],'fontsize',hDV.ftsz,'fontunits','normalized');
            grid(hDV.plotdata.pffot.axCBAR,'off');
            set(hDV.plotdata.pffot.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
            xlabel(hDV.plotdata.pffot.axCBAR,'Delta PP to slip [psi]')
            patchXs=linspace(0,1, hDV.cmapInfo.patchesInColorbar+1); % x coordinates of bounds of patches
            PatchYCoord= [0,1,1,0];% y coordinates of patches
            for j=1:hDV.cmapInfo.patchesInColorbar % cycle over each patch 
                PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
                hDV.plotdata.pffot.colorBarPatches(j,1)=patch(PatchXCoord,PatchYCoord,[.4 .4 .4],'parent',hDV.plotdata.pffot.axCBAR,'edgecolor','none') ;
            end
            
            
            %show stress regime in text
            hDV.plotdata.pffot.srtxt=uicontrol('parent',pl,'style','text','string','Stress Regime: Normal Faulting','position',[.55 .92 .4 .05],...
                'fontsize',14,'fontunits','normalized','backgroundcolor',[1 1 1]);
            
            %Mohr circle plot
            hDV.plotdata.pffot.ax2=axes('parent',pl,'position',[.55  .6  .4  .3],'color',hDV.colors.axsBgrnd,'XColor','k','YColor','k','fontsize',hDV.ftsz) ;
            hold(hDV.plotdata.pffot.ax2,'on');
            
            hDV.plotdata.pffot.mc1=plot(hDV.plotdata.pffot.ax2,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
            hDV.plotdata.pffot.mc2=plot(hDV.plotdata.pffot.ax2,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
            hDV.plotdata.pffot.mc3=plot(hDV.plotdata.pffot.ax2,0,0,'color',.5*[1 1 1],'linewidth',mohrCircleLinewidth);
            
                        
            % principal stress labels on mohr circle x axis
            hDV.plotdata.pffot.SigmaHhVLabel1=text([NaN],[NaN],{'\sigmaV'},'parent',hDV.plotdata.pffot.ax2,'color','c','HorizontalAlignment','center',...
                'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','bold','VerticalAlignment','baseline','clipping','on');
                        hDV.plotdata.pffot.SigmaHhVLabel2=text([NaN],[NaN],{'\sigmah'},'parent',hDV.plotdata.pffot.ax2,'color','c','HorizontalAlignment','center',...
                'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','bold','VerticalAlignment','baseline','clipping','on');
                        hDV.plotdata.pffot.SigmaHhVLabel3=text([NaN],[NaN],{'\sigmaH'},'parent',hDV.plotdata.pffot.ax2,'color','c','HorizontalAlignment','center',...
                'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','bold','VerticalAlignment','baseline','clipping','on');

            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.pffot.mflt(j)=plot(hDV.plotdata.pffot.ax2,0,0,'ko','markerfacecolor','k','markersize',hDV.plotdata.mohrCircleDotSize,'markeredgecolor','none');
            end
            
            hDV.plotdata.pffot.mfln =plot(0,0,'r','linewidth',2);
            axis(hDV.plotdata.pffot.ax2,'equal');
            xlabel(hDV.plotdata.pffot.ax2,'\sigma effective [psi]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pffot.ax2,'\tau [psi]','fontsize',hDV.ftsz) ;
            grid(hDV.plotdata.pffot.ax2,'on');
            
            
            %Stereonet plot
            hDV.plotdata.pffot.ax3=axes('parent',pl,'position',[.55  .08  .45  .45],'fontsize',hDV.ftsz) ;

            polarModified(0,0,'parent',hDV.plotdata.pffot.ax3) ;
            hold(hDV.plotdata.pffot.ax3,'on') ;
            for angleDipLines=15:15:75 % make line on s
                theta2=[0:pi/360:2*pi]';
                rho2=ones(size(theta2)).*angleDipLines./90;
                thisLine=polar(hDV.plotdata.pffot.ax3,theta2,rho2,':');
                set(thisLine,'color',[.5,.5,.5],'linewidth',0.5)
            end
            
            
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.snet(j) = polar(0,0,'parent',hDV.plotdata.pffot.ax3) ;
                hDV.plotdata.snetpoles(j) = polar(0,0,'o','parent',hDV.plotdata.pffot.ax3) ;
                set(hDV.plotdata.snet(j),'linewidth',2) ;
            end
            
            hDV.plotdata.snetcomp = surf(eye(2),eye(2),eye(2),'visible','off') ; % stereonet composite
            shading(hDV.plotdata.pffot.ax3,'interp');
            colormap(hDV.plotdata.pffot.ax3,flipud(hDV.cmapGYR)) ;
            set(hDV.plotdata.snetpoles(:),'visible','off','markersize',hDV.plotdata.mohrCircleDotSize) ;
            hDV.plotdata.sHmaxdir = polar(0,0,'parent',hDV.plotdata.pffot.ax3) ; %max hor stress direction
            set(findall(hDV.plotdata.pffot.ax3,'type','patch'),'visible','off')  % this prevents white patch objects from covering up stereonet poles if restacked
%              set(hDV.plotdata.pffot.ax3,'view',[0 -90]);
            
            uicontrol('parent',pl,'style','text','string','Stereonet Show:','position',[.55 .01 .15 .05],...
                'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1]);
            hDV.plotdata.pffot.popsnet = uicontrol('parent',pl,'style','popup','string','Fault Normals |Projected Curves |Normal Composite','position',[.7 .01 .2 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callsnetpop,hDV});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        case 'PROB. GEOMECH'
            %button to load stochastic data
            hDV.plotdata.pprob.hbutMCDat = uicontrol('parent',pl,'style','push','string','Load Distributions','position',[.05 .9 .19 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@calldatmc,hDV},'TooltipString','enter stochastic geomechanical uncertainties');
            %button to run analysis
            hDV.plotdata.pprob.hbutMC = uicontrol('parent',pl,'style','push','string','Run Analysis','position',[.26 .9 .19 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callrunmc,hDV},'backgroundcolor','r','foregroundcolor','w',...
                'TooltipString','run probabilistic geomechanical analysis');

            %plotting for colored risk cdf curves showing failure prob as
            %function of pressure
            hDV.plotdata.pprob.ax2=axes('parent',pl,'position',[.07  .3  .38  .45],'fontsize',hDV.ftsz,'color',hDV.colors.axsBgrnd) ;
            hold(hDV.plotdata.pprob.ax2,'on') ;
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.flinesprob(j)=plot(0,0,'linewidth',2,'color','k','parent',hDV.plotdata.pprob.ax2) ;
                hold(hDV.plotdata.pprob.ax2,'on') ;
            end
            xlim(hDV.plotdata.pprob.ax2,[0 100]) ;
            ylim(hDV.plotdata.pprob.ax2,[0 1]) ;
            grid(hDV.plotdata.pprob.ax2,'on');
            xlabel(hDV.plotdata.pprob.ax2,'\Delta Pore Pressure to Slip [psi]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pprob.ax2,'Probability of Fault Slip','fontsize',hDV.ftsz) ;
            
            %x axis control for risk curves plot
            uicontrol('parent',pl,'style','text','string','Max Delta PP [psi]:','position',[.17 .17 .2 .05],...
                'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1]);
            hDV.plotdata.pprob.htxtMaxDP=uicontrol('parent',pl,'style','edit','string','10000','position',[.37 .17 .1 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callmaxdpprob,hDV.plotdata.pprob.ax2,hDV});
            callmaxdpprob(hDV.plotdata.pprob.htxtMaxDP,[],hDV.plotdata.pprob.ax2,hDV);
            
            %plotting of sensitivity analysis bars
            hDV.plotdata.pprob.ax3=axes('parent',pl,'position',[.65  .1  .3  .35],'fontsize',hDV.ftsz) ;
            % Names of the Y axis ticks
            dat={'Vert Stress Grad' , 0 , 0 ;
                 'Shmin Gradient' , 0 , 0 ;
                 'SHmax Gradient' , 0 , 0 ;
                 'Pore Press  Grad' , 0 , 0 ;
                 'Strike of fault ' , 0 , 0 ;
                 'Dip of fault    ' , 0 , 0 ;
                 'SHmax Azimuth' , 0 , 0 ;
                 'Friction Coeff  ' , 0 , 0 ;
                };

            hDV.plotdata.pprob.names = dat(:,1) ;
            base_value=0; % The base value is where the y axis is centered
            low_vals  = cell2mat(dat(:,2)) ; 
            high_vals = cell2mat(dat(:,3)) ;
            [~ , ind]=sort(abs(low_vals-high_vals),'ascend');
            high_vals=high_vals(ind);
            low_vals=low_vals(ind);
            names= hDV.plotdata.pprob.names(ind);
            
            hDV.plotdata.pprob.barhigh = barh(high_vals,'r','parent','parent',hDV.plotdata.pprob.ax3);
            hold(hDV.plotdata.pprob.ax3,'on');
            hDV.plotdata.pprob.barlow = barh(low_vals,'r','parent','parent',hDV.plotdata.pprob.ax3);
            bh = get(hDV.plotdata.pprob.barhigh,'BaseLine');
            set(bh,'BaseValue',base_value);
            set(hDV.plotdata.pprob.ax3,'yticklabel',names);
            set(hDV.plotdata.pprob.ax3,'Ytick',1:length(names),'YTickLabel',1:length(names));
            set(hDV.plotdata.pprob.ax3,'yticklabel',names);
            set(gca,'XColor',.1*[1 1 1],'YColor',.1*[1 1 1]);
            
            grid(hDV.plotdata.pprob.ax3,'on');
            xlabel(hDV.plotdata.pprob.ax3,'\Delta Pore Pressure to Slip [psi]','fontsize',hDV.ftsz) ;
            title(hDV.plotdata.pprob.ax3,'Choose Fault for Sensitivity Analysis','fontsize',hDV.ftsz) ;
            
            %plotting of input variation display
            hDV.plotdata.pprob.ax4=axes('parent',pl,'position',[.65  .6  .3  .35],'fontsize',hDV.ftsz) ;
            % Names of the Y axis ticks

            hDV.plotdata.pprob2.names = dat(:,1) ;
            base_value=0; % The base value is where the y axis is centered
            
            hDV.plotdata.pprob2.barhigh = barh(high_vals,'r','parent','parent',hDV.plotdata.pprob.ax4);
            hold(hDV.plotdata.pprob.ax4,'on');
            hDV.plotdata.pprob2.barlow = barh(low_vals,'r','parent','parent',hDV.plotdata.pprob.ax4);
            bh = get(hDV.plotdata.pprob2.barhigh,'BaseLine');
            set(bh,'BaseValue',base_value);
            set(hDV.plotdata.pprob.ax4,'yticklabel',names);
            set(hDV.plotdata.pprob.ax4,'Ytick',1:length(names),'YTickLabel',1:length(names));
            set(hDV.plotdata.pprob.ax4,'yticklabel',names);
            set(gca,'XColor',.1*[1 1 1],'YColor',.1*[1 1 1]);

            grid(hDV.plotdata.pprob.ax4,'on');
            xlabel(hDV.plotdata.pprob.ax4,'Percent Deviation [%]','fontsize',hDV.ftsz) ;
            title(hDV.plotdata.pprob.ax4,'Variability in Inputs','fontsize',hDV.ftsz) ;
            
            %button to export CSV of CDF curves
            hDV.plotdata.pprob.hbutExportCDF = uicontrol('parent',pl,'style','push','string','Export CDF data','position',[.05 .05 .16 .05],...
                'fontsize',hDV.ftsz,'fontunits','normalized','callback',{@callExportCDF_CSV,hDV,hDV.plotdata.flinesprob},...
                'TooltipString','export CSV of Probability of Fault Slip given a pore pressure increase (CDF Curves)');
            hDV.plotdata.pprob.hbutShowInputHistograms = uicontrol('parent',pl,'style','push','string','Show Input Distributions','position',[.22 .05 .25 .05],...
                'fontsize',hDV.ftsz,'fontunits','normalized','callback',{@makeInputHistogramsProbabilistic,hDV},...
                'TooltipString','show blue histograms of input uncertainties');   % make figure of histograms button
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case  'PROB. HYDRO'
            
            %plotting for delta p histograms
            %function of pressure
            hDV.plotdata.hydprob.ax2=axes('parent',pl,'position',[.15  .3  .38  .45],'fontsize',hDV.ftsz,'color',hDV.colors.axsBgrnd) ;
            hold(hDV.plotdata.hydprob.ax2,'on') ;
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.dplinesprob(j)=plot(NaN,NaN,'linewidth',2,'color',hDV.colors.blue,'parent',hDV.plotdata.hydprob.ax2,'LineStyle','-') ;
                hDV.plotdata.greyCDFBackground(j)=plot(0,0,'linewidth',2,'color',hDV.colors.darkgrey,'parent',hDV.plotdata.hydprob.ax2,'LineStyle','--') ;
%                 hold(hDV.plotdata.hydprob.ax2,'on') ;
            end
            xlim(hDV.plotdata.hydprob.ax2,[0 100]) ;
            ylim(hDV.plotdata.hydprob.ax2,[0 1]) ;
            grid(hDV.plotdata.hydprob.ax2,'on');
            xlabel(hDV.plotdata.hydprob.ax2,'\Delta Pore Pressure on fault [psi]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.hydprob.ax2,{'Probability '},'fontsize',hDV.ftsz,'color','k') ;
            title(hDV.plotdata.hydprob.ax2,{'Probability of Exceedance';'of Pressure on fault'},'fontsize',hDV.ftsz,'color',hDV.colors.blue) ;
            
            %x axis control for pp curves plot
            uicontrol('parent',pl,'style','text','string','Max Delta PP [psi]:','position',[.17 .17 .2 .05],...
                'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1]);
            hDV.plotdata.hydprob.htxtMaxDP=uicontrol('parent',pl,'style','edit','string','10000','position',[.37 .17 .1 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callmaxdpprob,hDV.plotdata.hydprob.ax2,hDV});
            callmaxdpprob(hDV.plotdata.hydprob.htxtMaxDP,[],hDV.plotdata.hydprob.ax2,hDV);
            
            %button to load stochastic data
            hDV.plotdata.hydprob.hbutMCDat = uicontrol('parent',pl,'style','push','string','Load Distributions','position',[.05 .9 .19 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@calldatHydroMc,hDV});
            %button to run analysis
            hDV.plotdata.hydprob.hbutMC = uicontrol('parent',pl,'style','push','string','Run Analysis','position',[.26 .9 .19 .05],...
                'fontsize',12,'fontunits','normalized','callback',{@callRunHydro_mc,hDV},'backgroundcolor',hDV.colors.red,'foregroundcolor','w');
             %slider
            GW=0.005; HE = 0.06 ; WC = 0.18 ; WM = 1 - (WC + 3 * GW); PI = [.005 0.005];
            [hDV.hdsldr(3), hDV.hdsldr_txt(3)] = dateSldr(pl,PI,WM,HE,GW,hDV); %Slider #3, which belongs to probabilistic hydrology panel
            %button to export CSV of CDF curves
            hDV.plotdata.hydprob.hbutExportCDF = uicontrol('parent',pl,'style','push','string','Export Blue Curves','position',[.35 .005 .2167 .06],...
                'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callExportCDF_CSV,hDV,hDV.plotdata.dplinesprob},...
                'TooltipString','export CSV of Probability of Fault Pressure Exceedance');
            hDV.plotdata.hydprob.hbutShowInputHistograms = uicontrol('parent',pl,'style','push','string','Show Input Distributions','position',[0.572 .005 .205  .06],...
                'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callMakeInputHistogramsHydro,hDV},...
                'TooltipString','display histograms of hydrologic uncertainties');   % make figure of histograms button
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'INTEGRATED'
             % popup to choose what to show as labels on plot
            hDV.plotdata.pint.pop(1) = uicontrol('parent',pl,'style','popup','string',...
                'Choose Plot Labels|a) Fault Number|b) PP Change at fault [psi]|c) Fault Slip Potential',...
                'position',[.05 .85 .3 .1],'fontsize',14,'fontunits','normalized','callback',{@callintpop,hDV},'value',1);
            
%             %help button
%             hDV.plotdata.pint.buthelp = uicontrol('parent',pl,'style','push','string','Help',...
%                 'position',[.37 .9 .08 .06],'fontsize',14,'fontunits','normalized','callback',{@callinthelp,hDV});
    
            %main plot axes to map faults
            hDV.plotdata.pint.ax=axes('parent',pl,'position',[.07  .3  .38  .45],'fontsize',hDV.ftsz) ;
            
            %plotting of background surface
            Xmax=10; Ymax=Xmax; npts=100 ; %grid spacing for plotting
            [hDV.plotdata.pint.Xgrid , hDV.plotdata.pint.Ygrid] = meshgrid(linspace(0,Xmax,npts) , linspace(0,Ymax,npts)) ;
            hDV.plotdata.pint.Zgrid = 0*hDV.plotdata.pint.Xgrid ;
            hDV.plotdata.pint.hsf = surf(hDV.plotdata.pint.Xgrid , hDV.plotdata.pint.Ygrid , hDV.plotdata.pint.Zgrid , 'facealpha',.2,'parent',hDV.plotdata.pint.ax);
            hold(hDV.plotdata.pint.ax,'on') ; shading(hDV.plotdata.pint.ax,'interp');
            axis(hDV.plotdata.pint.ax,'equal'); axis(hDV.plotdata.pint.ax,'tight');
            view(hDV.plotdata.pint.ax,[0 90]) ;
            set(hDV.plotdata.pint.ax,'clim',get(hDV.plotdata.pflot.ax,'clim') ) ;
            %plotting colored faults placeholder
            for j=1:1:hDV.data.NFAULTSMAX
                hDV.plotdata.flinesintback(j)=plot(0,0,'linewidth',2,'color',hDV.colors.darkgrey,'parent',hDV.plotdata.pint.ax) ;
                hDV.plotdata.flinesint(j)=plot(0,0,'linewidth',3,'color','k','parent',hDV.plotdata.pint.ax) ;
                hDV.plotdata.flinesinttxt(j)=text(0,0,0,' ','fontsize',12,'fontunits','normalized','parent',hDV.plotdata.pint.ax,'clipping','on');
            end
            %set(hDV.plotdata.flinesinttxt(:),'visible','off') ;
            xlabel(hDV.plotdata.pint.ax,'x easting [km]','fontsize',hDV.ftsz) ;
            ylabel(hDV.plotdata.pint.ax,'y northing [km]','fontsize',hDV.ftsz) ;
            grid(hDV.plotdata.pint.ax,'on');
            
             % square for each well
            hDV.plotdata.pint.hwlocs=plot(NaN,NaN,'ks','linewidth',2,'markersize',8,'parent',hDV.plotdata.pint.ax) ;
            
            
            %plot fault pressure with time
            hDV.plotdata.pint.ax2=axes('parent',pl,'position',[.55  .52  .4  .35]) ;
            title(hDV.plotdata.pint.ax2,'Select Fault to Plot Pressures');
            hold(hDV.plotdata.pint.ax2,'on') ;            

            %time bars on the right plots
            hDV.plotdata.pint.timebar1=plot(hDV.plotdata.pint.ax2,NaN,NaN,'g--','linewidth',2,'displayname','Time Bar') ; 
            
            % make curves of pressure on fault
             for j=1:1:hDV.data.NFAULTSMAX
               hDV.plotdata.pint.PressureCurvesThruTime(j) = plot(0,0,'parent',hDV.plotdata.pint.ax2,...
               'linewidth',2,'HandleVisibility','off','color','b','marker','x');
             end
            
            %hDV.plotdata.pint.backar = area(0,0,'parent',hDV.plotdata.pint.ax2, 'FaceColor', [1 .7 .7] , 'EdgeColor' , [1 .7 .7],'displayname','Limit Prob. Variation');
            hDV.plotdata.pint.fltpressure = plot(hDV.plotdata.pint.ax2,0,0,'linewidth',2,'displayname','Deterministic Pressure at Fault') ;
            hDV.plotdata.pint.fltpressurelim = plot(hDV.plotdata.pint.ax2,0,0,'c','linewidth',2,'displayname','Deterministic PPF') ;
            hDV.plotdata.pint.fltpressurelimUp = plot(hDV.plotdata.pint.ax2,0,0,'--','linewidth',2,'displayname','50% P(slip)','color',[.7 .7 .7]) ;
            hDV.plotdata.pint.fltpressurelimDown = plot(hDV.plotdata.pint.ax2,0,0,'linewidth',2,'displayname','20% P(slip)','color',[.7 .7 .7]) ;
            hDV.plotdata.pint.fltpressurelimBottom = plot(hDV.plotdata.pint.ax2,0,0,':','linewidth',2,'displayname','1% P(slip)','color',[.7 .7 .7]) ;
            
            hDV.plotdata.pint.probHydroMaxMinPp = plot(hDV.plotdata.pint.ax2,NaN,NaN,'b+','visible','on','displayname','highest and lowest possible Pp','linewidth',3) ; % highest possible Pp
            hDV.plotdata.pint.probHydroMaxMinPpLabel=text([NaN;NaN],[NaN;NaN],[NaN;NaN],{' Highest Possible Pp ';' Lowest possible Pp '},'parent',hDV.plotdata.pint.ax2,'color','b',...
                'fontsize',hDV.ftsz,'fontunits','normalized','fontweight','normal','clipping','off','visible','on');

            xlabel('Time [years]','fontsize',hDV.ftsz) ; ylabel('Pressure Change at Fault Midpoint [psi]','fontsize',hDV.ftsz,'fontunits','normalized') ;
            grid(hDV.plotdata.pint.ax2,'on');% legend(hDV.plotdata.pint.ax2,'show','location','northoutside'); % would like to make fontunits normalized but don't know how to. 



            %controls for color bar range
            hDV.plotdata.pint.ppm=NaN; %placeholder for pressure margin numbers
            hDV.plotdata.pint.cmintxt = uicontrol('parent',pl,'style','edit','string','--','position',[.025 .12 .05 .03],...
                'fontsize',12,'fontunits','normalized','callback',{@callintcbartxt,hDV});
             hDV.plotdata.pint.cmaxtxt = uicontrol('parent',pl,'style','edit','string','--','position',[.425 .12 .05 .03],...
                'fontsize',12,'fontunits','normalized','callback',{@callintcbartxt,hDV});
            
            % plot colorbar axis for integrated
            hDV.plotdata.pint.axCBAR=axes('parent',pl,'position',[.05  .15  .4  .05],'fontsize',hDV.ftsz,'fontunits','normalized');
            grid(hDV.plotdata.pint.axCBAR,'off');
            set(hDV.plotdata.pint.axCBAR,'ytick',[],'XTickMode','auto', 'XTickLabelMode', 'auto');
            xlabel(hDV.plotdata.pint.axCBAR,'Fault Slip Potential')
            patchXs=linspace(0,1, hDV.cmapInfo.patchesInColorbar+1); % x coordinates of bounds of patches
            PatchYCoord= [0 1 1 0];% y coordinates of patches
            for j=1:hDV.cmapInfo.patchesInColorbar % cycle over each patch
                PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
                hDV.plotdata.pint.colorBarPatches(j,1)=patch(PatchXCoord,PatchYCoord,[.4 .4 .4],'parent',hDV.plotdata.pint.axCBAR,'edgecolor','none') ;
            end

           
           % these 2 lines allow me to turn off the FSP through time plot 
           visibility1='on';
           hDV.plotdata.pint.FSPThruTimeVisibility1=visibility1;
           
            %plot fsp with time on each fault
            hDV.plotdata.pint.ax4=axes('parent',pl,'position',[.55  .07  .4  .3],'visible',visibility1) ; % ,'fontsize',hDV.ftsz
            title(hDV.plotdata.pint.ax4,'FSP');
            hold(hDV.plotdata.pint.ax4,'on') ;            
            xlabel(hDV.plotdata.pint.ax4,'Time [years]','fontsize',hDV.ftsz) ; ylabel(hDV.plotdata.pint.ax4,'Fault Slip Potential','fontsize',hDV.ftsz) ;
            FSPYlims=[0,1];
            set(hDV.plotdata.pint.ax4,'ylim',FSPYlims);
            
            % background color patches
            hDV.plotdata.pint.FSPBackgroundPatches=zeros(length(hDV.cmapInfo.patchesInColorbar),1);
            PatchXs= [NaN,NaN,NaN,NaN];% coordinates of patches
            patchYs=PatchXCoord;
            for j=1:hDV.cmapInfo.patchesInColorbar % cycle over each patch
                hDV.plotdata.pint.FSPBackgroundPatches(j,1)=patch(PatchXs,patchYs,[0,0,0],'parent',hDV.plotdata.pint.ax4,'edgecolor','none','visible','off') ;
                
            end

            
%            % plot patches of Fault FSP through time
%             hDV.plotdata.pint.colorBarPatches=axes('parent',hDV.plotdata.pint.ax4,'fontsize',hDV.ftsz,'fontunits','normalized');
%             patchXs=linspace(0,1, hDV.cmapInfo.patchesInColorbar+1); % x coordinates of bounds of patches
%             PatchYCoord= [0 1 1 0];% y coordinates of patches
%             for j=1:hDV.cmapInfo.patchesInColorbar % cycle over each patch
%                 PatchXCoord = [patchXs(j),patchXs(j),patchXs(j+1),patchXs(j+1)] ;
%                 hDV.plotdata.pint.colorBarPatches(j,1)=patch(PatchXCoord,PatchYCoord,[.4 .4 .4],'parent',hDV.plotdata.pint.axCBAR,'edgecolor','none') ;
%             end
            
            faultTextLabelsMax=9;
            % make curves of FSP vs time for each fault
             for j=1:1:hDV.data.NFAULTSMAX
               hDV.plotdata.pint.FaultFSPCurvesThruTime(j) = plot(0,0,'parent',hDV.plotdata.pint.ax4,... % curve
               'linewidth',2,'HandleVisibility','off','color','k','visible',visibility1,'marker','x');
               hDV.plotdata.pint.FSPtextLabelTime2(1:faultTextLabelsMax,j)= text(nan(faultTextLabelsMax,1),nan(faultTextLabelsMax,1),num2str(j),...
                   'verticalalignment','bottom','color','k','parent',hDV.plotdata.pint.ax4,'clipping','on','visible',visibility1); % fault number text
             end
               %time bars on the right plots
            hDV.plotdata.pint.timebar2=plot(hDV.plotdata.pint.ax4,NaN,NaN,'g--','linewidth',2,'displayname','Time Bar FSP thruTime','visible',visibility1) ; 
             grid(hDV.plotdata.pint.ax4,'on');
            %button to export CSV of FSP vs time curves
            hDV.plotdata.pprob.hbutExportFSPData = uicontrol('parent',pl,'style','push','string',{'Export'},'position',[.91  .37  .09  .05],...
                'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callExportFSP_CSV,hDV},'visible',visibility1);  
             
           % summarize Fault FSP through time button?
           hDV.plotdata.pprob.summaryPlotButton = uicontrol('parent',pl,'style','push','string',{'Summary Plots'},'position',[.2  .8 .19 .05],...
                'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callRunSummaryPlot,hDV},'backgroundcolor',hDV.colors.red,'foregroundcolor','w',...
                'TooltipString','(This is Slow!!) Click to generate a summary plot of FSP and pressure through time');
            
            %slider
            GW=0.005; HE = 0.06 ; WC = 0.18 ; WM = 1 - (WC + 3 * GW); PI = [.005 0.005];
            [hDV.hdsldr(2), hDV.hdsldr_txt(2)] = dateSldr(pl,PI,WM,HE,GW,hDV); %Slider #2 belongs to INTEGRATED panel
            
             %button to export CSV of Pressure vs time curves
            hDV.plotdata.pprob.hbutExportFaultPres = uicontrol('parent',pl,'style','push','string',{'Export'},'position',[.91  .95  .09  .05],...
                'fontsize',hDV.ftsz-1,'fontunits','normalized','callback',{@callExportPP_CSV,hDV});


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%

    function callmaxcpf(src,~,hDV)
        hax  = hDV.plotdata.pflot.ax;
        hax2 = hDV.plotdata.pflot.ax2;
        hax3 = hDV.plotdata.pint.ax;
        
        val = str2double(get(src,'string'));
        cmin = get(hax,'clim'); cmin=cmin(1);
        if val<=0 || isnan(val) || val<cmin
            val=cmin+100;
        end
        set([hax ; hax3],'clim',[cmin val]);
        ylim(hax2,[0 val]) ;
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%


    function callchkmohr(src,~,hDV)
        if get(src,'Value')
            set(hDV.plotdata.pflot.mc1(:),'visible','on');
            set(hDV.plotdata.pflot.mc2(:),'visible','on');
            set(hDV.plotdata.pflot.mc3(:),'visible','on');
        else
            set(hDV.plotdata.pflot.mc1(:),'visible','off');
            set(hDV.plotdata.pflot.mc2(:),'visible','off');
            set(hDV.plotdata.pflot.mc3(:),'visible','off');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%
    function callsnetpop(src,~,hDV)
        switch get(src,'value')
            case 1 % selected faults
                set(hDV.plotdata.snet(:),'visible','off');
                set(hDV.plotdata.snetpoles(:),'visible','on');
                set(hDV.plotdata.snetcomp,'visible','off');
                callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV) % plot faults on stereonet that are selected
            case 2 % selected faults
                set(hDV.plotdata.snet(:),'visible','on');
                set(hDV.plotdata.snetpoles(:),'visible','off');
                set(hDV.plotdata.snetcomp,'visible','off');
                callbackFaultsSelected(hDV.plotdata.ListboxFaultSelector,[],hDV) % plot faults on stereonet that are selected
            case 3 % all orientations
                set(hDV.plotdata.snet(:),'visible','off');
                set(hDV.plotdata.snetpoles(:),'visible','off');
                set(hDV.plotdata.snetcomp,'visible','on');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%






end