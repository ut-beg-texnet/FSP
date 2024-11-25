% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface
%
% **** Completely rewritten by Suvrat Lele for FSP 2.0 April 2018 ****
% Subsurface Mechanics and Induced Seismicity, Drilling & Subsurface
% ExxonMobil URC


% Mohr's Circle Calculation 3D
%
%Inputs:
% Cell array indatacell:
%     1 Sig0 = Vector of initial total stresses (assume row vector) 
%     2 ignoreParameter4 (Was ixSv, the index of vertical stress)
%     3 p0   = Initial Reference pore pressure corresponding to Sig0
%     4 strikes of faults in degrees
%     5 dips of faults in degrees
%     6 SHdir= Maximum horizontal stress direction
%     7 dp   = Pressure pertubation at each fault
%     8 mu   = friction coefficients for faults
%     9 biot coeff
%     10 nu Poisson's ratio
%     11 aphi if 11th cell exists
% hDV  = main data structure
%
% Outputs:
% failout    = Delta pore pressure to failure; only output used in Monte Carlo
% outs.ppfail = Delta pore pressure to failure for fault (includes Poisson's effect)
% outs.cff    = Couloumb failure function (effective shear stress to fail - not changed for Poisson's effect)
% outs.scu    = Shear capacity utilization
% C1, C2, C3 = Plotted circles
% sig_fault  = Effective normal stress projected onto fault
% tau_fault  = Effective shear stress projected onto fault


%%%%% IMPORTANT! Order of inputs must match ppfail_MC and other MC calcs
function [failout,outs,C1,C2,C3,sig_fault,tau_fault]=mohrs_3D(indatacell,hDV)

Sig0 = indatacell{1} ;
p0 = indatacell{3} ;
strike = indatacell{4} ;
dip = indatacell{5} ;
SHdir = indatacell{6} ;
dp = indatacell{7} ;
mu = indatacell{8} ;
biot = indatacell{9} ;
nu = indatacell{10} ;

% If 11th input exists, calculate stresses using appropriate A-phi model
%
if numel(indatacell)>10
    % Use fault mu as ref mu in FSP 2.0
    % Deterministic tab may analyze multiple faults at the same time, but all
    % faults have the same mu, so use the first one, mu(1)
    % Monte Carlo runs only one fault at one time with only one set of
    % stress and fault parameters, so mu(1) is the only value in that case
    %
    % Either mu or shmin is used accoring to aphi.use value indicating whether
    % Modified A-Phi model is to be used
    %
    APhi = indatacell{11};
    [SH,Sh] = getHorFromAPhi(APhi,mu(1),Sig0(2),Sig0(1),p0,hDV.data.stress.aphi.use);
    Sig0(2) = Sh;
    Sig0(3) = SH;
end


if any(isnan(dp));disp(['found DP''s =NaN, heads up, setting =0']);end
dp(isnan(dp))=0; % this was causing an error on some falults. I FRW am not certain that this is the correct fix but it seems to stop errors being thrown. 



% Calculate delta pore pressure to slip (failure) and the resolved shear and 
% normal stresses on fault - Poisson's ratio effect is included in the function below
%
% All faults have same mu in FSP 2.0, so use the first one
%
[failout, tau_fault, sig_fault] = calc_ppfail(Sig0, SHdir, p0, biot, nu, mu(1), dp, strike, dip);



% The following calculates various output variables that are only used in 
% deterministic analysis
% Monte Carlo only uses the first output, failout above, so do not run code 
% below if the extra outputs are not requested when function is called
%
if nargout > 1
    
    outs.ppfail = failout;
    outs.cff    = tau_fault-mu.*sig_fault;
    outs.scu    = tau_fault./(mu.*sig_fault);

    % The code below for Mohr circle plotting is unchanged from FSP 1.0,
    % except moving ixSv calculation here (ixSv removed as input above)
    %
    %compute poroelastic total stresses
    %
    N=length(strike) ; 
    Sig0sorted = sort(Sig0,'descend');
    ixSv = find(Sig0sorted==Sig0(1)); %index of vertical stress - moved here by SPL
    ixSv = ixSv(1); %use only one of the indices if there is more than one
    Sig = kron(Sig0sorted,ones(N,1));
    Ds  = kron(biot*(1-2*nu)/(1-nu)*dp,[1 1 1]);
    Ds(:,ixSv) = 0 ;
    Sig = Sig + Ds ;

    a = linspace(0,pi); c=exp(1i*a);
    C1 = zeros(N,length(a)) ; 
    C2 = zeros(N,length(a)) ; 
    C3 = zeros(N,length(a)) ; 

    R1= .5*(Sig(:,1)-Sig(:,3));
    R2= .5*(Sig(:,2)-Sig(:,3));
    R3= .5*(Sig(:,1)-Sig(:,2));

    for k=1:1:N
        C1(k,:)=R1(k)*c+(Sig(k,1)+Sig(k,3))/2-(p0+dp(k)) ;
        C2(k,:)=R2(k)*c+(Sig(k,2)+Sig(k,3))/2-(p0+dp(k)) ;
        C3(k,:)=R3(k)*c+(Sig(k,1)+Sig(k,2))/2-(p0+dp(k)) ;
    end
    %
    % end code for Mohr circle plotting
    
end % nargout > 1



    % Function to calculate delta pore pressure to slip (failure)
    % Supports analysis for multiple faults at the same time
    % Vectorized code: dp, str, dip are vectors with one value per fault
    % Friction coeff mu is same for all faults and is a scalar
    %
    % Poisson's ratio effect is included
    %
    function [ppfail, tau_fault, sig_fault] = calc_ppfail(Sig0, az, p0, biot, nu, mu, dp, str, dip)

    % Cos and sin of azimuth for use below
    %
    cos_az = cos(az*3.1415926d0/180.d0);
    sin_az = sin(az*3.1415926d0/180.d0);

    % Cos and sin of strike and dip angles for use below
    % Note that str and dip are arrays for multiple faults - vectorized code
    %
    cs = cos(str*3.1415926d0/180.d0);
    ss = sin(str*3.1415926d0/180.d0);
    cd = cos(dip*3.1415926d0/180.d0);
    sd = sin(dip*3.1415926d0/180.d0);
    
    % Total stresses from input data; new variables created for readability of code
    % 
    Svert = Sig0(1);
    shmin = Sig0(2);
    sHmax = Sig0(3);
    
    
    % Factor based on nu for horizontal stress change
    % Biot is currently hardcoded to 1 (outside this function), and
    % probably good for fault stress analysis, similar to Abaqus where
    % effective stress is always total-pp even when Biot effect is included
    % via *Porous bulk moduli option
    %
    f = biot*nu/(1-nu);  
    
    
    % Effective stresses at current pressure, compressive positive notation
    % dp is also an array, which makes s11, s22, s33, s12 also arrays
    %
    % Note that factor f is used for dp because dp is pressure change due to
    % injection upto present time - horizontal stresses will have Poisson's
    % ratio effect (and Biot coeff, if other than 1; but see comments above)
    % However, p0 is initial pressure before injection, so factor f is not needed
    %
    s11 = shmin*cos_az*cos_az + sHmax*sin_az*sin_az - p0 - f*dp;
    s22 = shmin*sin_az*sin_az + sHmax*cos_az*cos_az - p0 - f*dp;
    s33 = Svert - p0 - dp;
    s12 = (sHmax - shmin)*cos_az*sin_az;
    
   
    % Components of unit normal vector to fault planes (multiple faults analyzed together)
    %
    n1 = sd .* cs;
    n2 = -sd .* ss;
    n3 = cd; % However, note that n3 is not used below, since by def n3 = 1-sqrt(n1^2+n2^2)


    
    % Shear and normal stresses resolved on fault
    % Note that tau_fault is absolute magnitude, but sig_fault is signed
    % value (which should be positive for most cases, but, with high pore
    % pressure, resolved normal stress on fault can become negative any of 
    % the principal stresses is negative)
    %
    tau_fault = sqrt((n2.^2 .*(s12.^2 -(-1+n2.^2).*(s22-s33).^2) ...
                  - n1.^4 .*(s11-s33).^2  ...
                  + 4*n1.^3 .*n2 .*s12.*(-s11+s33) ... 
                  + 2*n1 .* n2 .*s12 .*(s11 +s22 - 2*n2.^2 .*s22 + 2*(-1+n2.^2) .*s33) ...
                  + n1.^2 .*(s11.^2 +(1-4*n2.^2).*s12.^2 - 2*s11.*(n2.^2 .*(s22-s33) +s33) +s33.*(2*n2.^2 .*(s22-s33) +s33))));
              
    sig_fault =  2*n1.*n2.*s12 + n1.^2 .*(s11-s33) + n2.^2 .*(s22-s33) +s33;  % signed sig_fault

    % Mobilized friction coeff with present stress and pressure, used below
    % to determine if pressure to slip should be positive or negative
    %
    mobmu = tau_fault./sig_fault;
    mobmu(isnan(mobmu)) = 99.99;  % a high number, if sig_fault is zero and mob_mu is NaN
    mobmu(mobmu<0)      = 99.99;  % again, set to high number if sig_fault is negative and so mobmu is negative

    
    % Coefficients of quadratic equation A*dp^2 + B*dp + C = 0
    % such that the mobilized friction coefficient on fault = mu
    % i.e. dp that will bring the fault stresses to failure point
    %
    C = -4*(1+mu^2)*n1.^3 .* n2.*s12.*(s11-s33) ...
        - (1+mu^2)*n1.^4 .* (s11-s33).^2 ...
        - (1+mu^2)*n2.^4 .* (s22-s33).^2 ...
        - mu^2 * s33.^2 ...
        + 2*n1 .* n2.*s12 .* (s11 + (1-2*(1+mu^2)*n2.^2).*s22 + 2*(1+mu^2)*(-1+n2.^2).*s33) ...
        + n2.^2 .* (s12.^2 + s22.^2 - 2*(1+mu^2)*s22.*s33 + (1+2*mu^2)*s33.^2) ...
        + n1.^2 .* (s11.^2 + (1-4*(1+mu^2)*n2.^2).*s12^2 - 2*(1+mu^2)*s11.*(n2.^2 .*(s22-s33) +s33) + s33.*(2*(1+mu^2)*n2.^2 .*(s22-s33) +s33 +2*mu^2 .*s33));

    B = 2*( 2*(-1+f)*(1+mu^2)*n1.^3 .* n2.*s12 ...
            + 2*n1.*n2 .* (-(1+mu^2)*(-1+n2.^2) + f*(-1+(1+mu^2)*n2.^2)).*s12 ...
            + (-1+f)*(1+mu^2)*n1.^4 .*(s11-s33) ...
            + (-1+f)*(1+mu^2)*n2.^4 .*(s22-s33) ...
            + mu^2 *s33 ...
            + n2.^2 .*((1-f+mu^2)*s22 + (-1 +f -2*mu^2 +f*mu^2)*s33) ...
            + n1.^2 .* ((-(1+mu^2)*(-1+n2.^2) + f*(-1+(1+mu^2)*n2.^2)).*s11 + (-1 +f)*(1+mu^2)*n2.^2 .*(s22-2*s33) +(-1+f-2*mu^2 +f*mu^2).*s33) );

    A = -mu^2 *(1+(-1+f)*n1.^2 +(-1+f)*n2.^2).^2 ...
        - (-1+f)^2 *(n1.^4 + n2.^2 .*(-1+n2.^2) + n1.^2 .*(-1+2*n2.^2));

    Bsq_minus_4AC = B.*B -4*A.*C;

    ppfail1 = (-B -sqrt(Bsq_minus_4AC))./(2*A);  % Two solutions of above
    ppfail2 = (-B +sqrt(Bsq_minus_4AC))./(2*A);  % quadratic equation
    
    
    % For some cases, with low Poisson's ratio, there may be no solution
    % to above equation, indicated by Bsq_minus_4AC <0, since its 
    % sqrt will be imaginary
    %
    % for such cases set ppfail to the horizontal distance from failure line 
    % (i.e. value for nu=0.5)
    % Assign this to ppfail1 and ppfail2 with opposite signs; mobmu
    % relative to mu below will determine correct sign (since there is no
    % solution this is just a placeholder anyway)
    %
    % Cannot set all these to same value such as zero becasue probability
    % distribution calculation (outside this function) has some issues
    %
    no_solution = Bsq_minus_4AC <0;  % logical array
    ppfail_horiz_dist = sig_fault - tau_fault/mu;
    ppfail1(no_solution) = -ppfail_horiz_dist(no_solution); 
    ppfail2(no_solution) =  ppfail_horiz_dist(no_solution); 

    
    
    ppfail = nan(size(dip));  %Initialize ppfail to NaN
    
    
    % Loop over faults (this will be very difficult to vectorize, and to read/debug if vectorized)
    %
    for k=1:length(dip)
        
        if mobmu(k) < mu   % If fault is "below failure line"
            
            if ppfail1(k)>0 && ppfail2(k)>0
                % If both roots are positive then choose smaller one
                if ppfail1(k)<ppfail2(k)  ppfail(k)=ppfail1(k);  else ppfail(k)=ppfail2(k);  end 
                
            elseif ppfail1(k)<0 && ppfail2(k)<0
                % Both roots negative - probably indicates some error
                % There has to be one positive solution, since at some
                % pressure all Mohr circles will move enough to the left 
                % 
                edlgBox=errordlg('Pressure to slip calculation error','Calculation Error');
                centerFigure(hDV.hfig,edlgBox);
                
            else
                % Now one should be negative and the other positive, choose the positive one
                if ppfail1(k)>0  ppfail(k)=ppfail1(k);  else ppfail(k)=ppfail2(k);  end
            end
            
        elseif mobmu(k) > mu  % fault is "above failure line"
            
            if ppfail1(k)<0 && ppfail2(k)<0
                % If both roots are negative then choose the one with smaller *magnitude*
                if ppfail1(k)<ppfail2(k)  ppfail(k)=ppfail2(k);  else ppfail(k)=ppfail1(k);  end
                
            elseif ppfail1(k)>0 && ppfail2(k)>0
                % Both roots positive, even though fault initially above failure line
                % This happens in a few cases with low nu due to nonlinear stress path 
                % In this case, output horizontal distance
                %
                ppfail(k) = sig_fault(k) - tau_fault(k)/mu;
                
            else
                % Now one should be negative and the other positive, choose the negative one
                if ppfail1(k)>0  ppfail(k)=ppfail2(k);  else ppfail(k)=ppfail1(k);  end
            end
            
        else % mobmu == mu (very unlikely due to finite numerical precision)
            ppfail(k) = 0;
        end
        
    end % loop over faults

    end % function deltaPP_slip


end % function mohr_3D