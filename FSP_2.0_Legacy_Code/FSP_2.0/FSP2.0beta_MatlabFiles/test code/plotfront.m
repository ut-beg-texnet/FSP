%%% start reset memory
clc ;clear all ;close all ; close force ;

rho=1000 ; g=9.8 ; 

%%% define basic geometry dimensions
Xmax=15; Ymax=Xmax;
npts=200 ; %grid spacing for plotting
xwell=10; ywell=xwell; %well location

[Xgrid , Ygrid] = meshgrid(linspace(0,Xmax,npts) , linspace(0,Ymax,npts)) ;

difflocs=Xgrid+1i*Ygrid - (xwell+1i*ywell) ;
R=abs(difflocs) ; %compute radial distances

S=5.4e-5; %m^2/sec storativity (dimensionless) 
T=4.5e-6; %m^2/sec transmissivity

Q=2*6.7e6/30/24/3600*1e-3; %m^3 per sec injection rate
t=10 ; %years, injection period

%%% Presuure calcultion
P=pfront(R*1e3,t*365*24*3600,Q,S,T,rho,g) ;
P = P+5*randn(size(P)) ; %add some noise

%%% Plot pressure field
figure; surf(Xgrid , Ygrid , P , 'facealpha',.8) ; shading interp ;
view([0 90]) ; xlabel('x easting [km]') ; ylabel('y northing [km]') ;  hold on ;
plot3(xwell,ywell,max(P(:)),'kx','linewidth',2,'markersize',10) ; colorbar ;
zlabel('Pressure front [bars]') ; axis equal ; axis tight;
caxis([5 40]) ;


%%% Faults definitions
nfaults = 100 ; %number of faults
xf=Xmax*rand(nfaults,1);
yf=Ymax*rand(nfaults,1);
thf=180*rand(nfaults,1); %fault orientation

% calculate pore pressure for each fault from the grid
dx=Xmax/npts ; dy=Ymax/npts ;
ppf=0*xf;
for k=1:1:nfaults
    ppf(k)=P(ceil(xf(k)/dx),ceil(yf(k)/dy));
end

% define fault locations/geometry
lenf=1 ;  
fstart = xf+1i*yf - lenf/2*exp(1i*thf*pi/180) ;
fend   = xf+1i*yf + lenf/2*exp(1i*thf*pi/180) ;

% calculate strike alpha
alphf = thf ;
alphf(thf>90)=180-thf(thf>90); 

%%% insitu stress and friction coefficient for each fault
sig=[250 200 100]; muf=0*thf+0.6;

%%% compute fault failure margin for each fault and plot Mohr's circle
figure; ax=axes;
cfc=mohrs_simple(sig,alphf,ppf,muf,ax) ;

%define a colormap
cmap = autumn(nfaults) ; 

% plot color coded faults
figure; hold on ;
for k=1:1:nfaults
    h=plot(real([fstart(k) fend(k)]) , imag([fstart(k) fend(k)]) ) ;
    cv=cfc;
    fr = (cv(k)-min(cv))/(max(cv)-min(cv)) ; 
    ix= 1+floor(fr*(nfaults-1));
    set(h, 'color', cmap(ix,:));
end
axis tight; axis equal;  grid on ; %whitebg([0 0 0])
axis([0 Xmax 0 Ymax]) ; 

% plot stress margin as a function of strike
figure; 
plot(alphf , cfc,'or') ; hold on ; xlabel('\alpha [deg]') ; ylabel('\sigma magin CFC [bars]') ;   
grid on;

% failed faults fraction and include on plot
f_failed = sum(cfc<0)/length(cfc) ; 
text(45,0.8*max(cfc),['\bf Fraction of failed faults: ',num2str(100*f_failed),'%'],'fontsize',12); 




