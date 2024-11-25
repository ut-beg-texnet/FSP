% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

% modified by Rall Walsh to take variable injection rates
% differences: t is now the year of interest (ie 2019), not the duration of
% injection in years or seconds (ie 5)
% Q is now the data structure known as: thisWellDatenumBarrelsPerDay for
% the given well {access with curly brackets}
% not one constant rate Q.
% with help and code from Matt Weingarten

%Radially Symmetric finite thickness pressure front calculation
%Storativity (S) and Transmissivity (T) calculated externally.
% T  = transmissivity read in, assume it is the same in both X and Y directions.
% rho = fluid density in SI (1000), g = accel due to gravity in SI (9.81)
% r= distance of calculation from well in meters.
% yearsToCalculate is in decimal years

function delta_p_Bars=pfront(r,yearsToCalculate,thisWellDatenumBarrelsPerDay,S,T,rho,g)
if nargin==2
    S=5.4e-5; %m^2/sec storativity
    T=4.5e-6; %m^2/sec transmissivity
    % Q=6.7e6/30/24/3600*1e-3; %m^3 per sec injection rate - the function now takes barrels per day and converts to Q below
end

% throw warning if not an integer
if (abs(round(yearsToCalculate)-yearsToCalculate)) > eps('double')
%     .. not integer... throw warning
disp(' right now pfront is set up to take integer years (ie. 2012), and this isn''t one')
end

% convert r from row vector or matrix to column vector if necessary -
% switches back at the end
r_reshaped=false;
[rowr,colr]=size(r);
if colr~=1 % is r a column vector?
    r_reshaped=true;
    r=reshape(r,[],1); % make into  1 column vector, will reshape it back at the end
end

%loop over time
for yearCycleCount=1:length(yearsToCalculate);% cycle over yearOfInterest here
    
    yearOfInterest=yearsToCalculate(yearCycleCount); % year of slider bar/inquiry
    
    thisWellDatenumBarrelsPerDayUpTilNow=thisWellDatenumBarrelsPerDay(thisWellDatenumBarrelsPerDay(:,1)<=datenum(yearOfInterest,1,1,0,0,0),:);  % crop data to everything before this year
    
    if size(thisWellDatenumBarrelsPerDayUpTilNow,1)<1 % is the year in question before injection started in this well?
        if r_reshaped
            delta_p_Bars=zeros(rowr,colr);
        else
            delta_p_Bars(1,yearCycleCount)=0;
        end
        continue  % then don't bother with this well-year combination, because the year in question is before injection started in this well, return zeros
    end
    
    % put last datapoint at the beginning of year that slider bar is at
    if     size(thisWellDatenumBarrelsPerDayUpTilNow,1) < size(thisWellDatenumBarrelsPerDay,1)  % then there were datapoints after time of inquiry that were cropped.
        thisWellDatenumBarrelsPerDayUpTilNow(end+1,1:2)=[datenum(yearOfInterest,1,1,0,0,0),thisWellDatenumBarrelsPerDay(size(thisWellDatenumBarrelsPerDayUpTilNow,1)+1,2)];
    else % there was not a data point at the beginning of the year where the slider is, so add a zero datapoint.
        thisWellDatenumBarrelsPerDayUpTilNow(end+1,1:2)=[datenum(yearOfInterest,1,1,0,0,1),0];
    end
    
    tdays = thisWellDatenumBarrelsPerDayUpTilNow(:,1)-thisWellDatenumBarrelsPerDay(1,1); % code requires first time datapoint be zero time, so subtract first datapoint time from all to get elapsed time.
    t = tdays .*24 .*3600 ; % elapsed time in days to time in seconds.
    Q = thisWellDatenumBarrelsPerDayUpTilNow(:,2) .* 1.84013e-6 ;% 1 oil barrel per day = 1.84013e-6 cubic meters/second according to google
    
    % Initialize Well Function Parameters
    u = zeros(length(r),1);  % u matrix: input equation into well-function [Dimensionless]
    well_func = u;     % Well-function matrix: From Allen (1954) and Hastings (1955) where u varies greatly [Found in Practical Approximations of Well Function (Srivastava and Guzman-Guzman)]
    tstep_sum = u;     % Place-holder matrix for each timestep to be summed
    
    for ii=2:length(t) % cycle over times (in seconds)
        u(:,1) = ((((r(:,1)).^2 ).*T).*S)./(4.*T.*T.*(max(t)-t(ii-1,1))); % u matrix calculation at each x,y grid location...
        well_func(:,1) = expint(u(:,1)); % well function
        tstep_sum = tstep_sum + well_func.*(Q(ii,1)-Q(ii-1,1)); % sum previous tstep_sum with new well_function and injection rate step. Note: If Q has not changed, tstep_sum unchanged.
    end
    head = tstep_sum.*(1/(4*pi*T)); % Final head distribution @ r,time
    delta_p_Pascals = (head.*rho.*g) ;  % in pascals
    
    if r_reshaped % reshape back to shape of r (if necessary)
        delta_p_Pascals =  reshape(delta_p_Pascals,rowr,colr);
        delta_p_Bars=delta_p_Pascals.*1E-5; % pascals to Bars , and make row vector
    else % r wasn't reshaped, maybe multiple timesteps
        delta_p_Bars(:,yearCycleCount)=delta_p_Pascals.*1E-5; % pascals to Bars , and make row vector
    end
    
    %    delta_p_PSI=delta_p_Pascals ./6894.76;  %change it to PSI:  1 PSI = 6894.76 Pascals
    delta_p_Bars(isnan(delta_p_Bars))=0;
    delta_p_Bars(isinf(delta_p_Bars))=0;
end
