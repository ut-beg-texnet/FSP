function p=pfront(r,t,Q,S,T,rho,g)

%assuming SI inputs

if t<=0
    p=0*r;
else
    
    if nargin==2
        S=5.4e-5; %m^2/sec storativity
        T=4.5e-6; %m^2/sec transmissivity
        Q=6.7e6/30/24/3600*1e-3; %m^3 per sec injection rate
    end
    
    A1 = rho*g*Q/4/pi/T ;
    u = r.^2*S/(4*T*t) ;
    
    p = A1*expint(u)*10^-5 ; %bars
    
end
end