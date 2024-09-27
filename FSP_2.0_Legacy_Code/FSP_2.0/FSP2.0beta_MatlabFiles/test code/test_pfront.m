r=linspace(.1,20) ; t=[5 10 15 20] ;

S=5.4e-5; %m^2/sec storativity
T=4.5e-6; %m^2/sec transmissivity
Q=6.7e6/30/24/3600*1e-3; %m^3 per sec injection rate
  rho = 1000;  g = 9.81  ;

figure ; hold all ;


for k=1:1:length(t)
    p=pfrontOldSaved(r*1e3,t(k)*365*24*3600,Q,S,T,rho,g) ;
    plot(r,p) ;
end
 
 
 
 xlabel('Distance from well [km]') ; ylabel('Pressure [bars]') ; grid on;
 ylim([0 30]) ; 
 legend(num2str(t')) ; 