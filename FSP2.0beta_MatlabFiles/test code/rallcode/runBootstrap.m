

clearvars 
close all

load('BootstrapParametersForDarren.mat')
% save('BootstrapParametersForDarren.mat','S1Absolutes','S2Absolutes','S3Absolutes','Pps','manyMus','a445','b445','c445')
% these numbers are at 5 km depth
% data in MPA
% each above is a 10,000 row by 1 column vector
% a,b,and c are described as such: 
% I think there is an issue with this: lower hemisphere vr right hand rule,
% but I'm still wrapping my head around it. 
% a=trend of S1...exception when S1 is vertical a=trend SHmax-90 degrees
% b=-plunge of S1 (plunge is angle down from horizontal)
% c=rake of S2, 0 if S1 or S3 is vertical, 90 if S2 is vertical
% the 445 in a445 just means it's a vector of 10,000 a vaules, etc. 

% this version has no uncertainty in fault orientation. 
% to add that, you'd simply vary the fault orientation going into each of
% the NBOOT iterations. 


NBOOT=10000; % number of bootstraps
faults=[23,12;140,85;160,80]; % strike,dip, each as columns
% dip is 90 degrees to the right of strike
fontSize=16;

faultIdx390=1:1:size(faults,1); % indexes of faults (counting)
nfaults=length(faultIdx390); % number of faults
strikes=ones(NBOOT,nfaults);% preallocate
dips=ones(NBOOT,nfaults);% preallocate
failureCriteria=1; % 1 means PPF, 3 means CFF, 2 will be SCU but hasn't been calculated yet
distanceToFailureResult=zeros(NBOOT,nfaults) ; % preallocate
cdfLineDataThisAreaYData=zeros(NBOOT+1,nfaults); % preallocate
cdfLineDataThisAreaXData=zeros(NBOOT+1,nfaults); % preallocate

for faultCountIndex=faultIdx390  % cycle over each fault
    for bootstrapCountIndex=1:NBOOT % cycle over bootstraps 
% this calls the function stress_on_faults9RunNBOOTTimes 
     [distanceToFailureResult(bootstrapCountIndex,faultCountIndex)]=-stress_on_faults9RunNBOOTTimes(S1Absolutes(bootstrapCountIndex,1),...
         S2Absolutes(bootstrapCountIndex,1),S3Absolutes(bootstrapCountIndex,1),Pps(bootstrapCountIndex,1),manyMus(bootstrapCountIndex,1),...
         a445(bootstrapCountIndex,1),b445(bootstrapCountIndex,1),c445(bootstrapCountIndex,1),faults(faultCountIndex,1:2),failureCriteria) ;
  
    end
    % this line calculates the empirical cumulative distribution function
    % from the distanceToFailureResult for each fault
     [cdfLineDataThisAreaYData(:,faultCountIndex),cdfLineDataThisAreaXData(:,faultCountIndex)] = ecdf(distanceToFailureResult(:,faultCountIndex));  % this line calculates the cdf curve coordinates

end




figure % plot data as PDF
hold on

    % simple histogram of results
% hist(distanceToFailureResult )%
EDGES=0:2:max(max(distanceToFailureResult));
N = histc(distanceToFailureResult,EDGES); % calcuylate histogram bin heights
plot(EDGES,N./NBOOT) % plot edges, binheights divided by NBOOT (should be bin centers, being lazy)


title('PDF of results for each fault','fontSize',fontSize)
xlabel('MPa','fontSize',fontSize)
ylabel('probability','fontSize',fontSize)
set(gcf,'color','w')
set(gca,'fontSize',fontSize)


 figure % plot CDF lines
 hold on
 % plot each CDF curve
plot(cdfLineDataThisAreaXData,cdfLineDataThisAreaYData,'k')
ylabel('probability','fontSize',fontSize)
xlabel('pore pressure to whatever failure criteria','fontSize',fontSize)
set(gca,'fontSize',fontSize)
title('cdf curve','fontSize',fontSize)
set(gcf,'color','w')
% slice curves at X percent, and plot red X there
probabilityToQueryAt=0.1; % 10 %
for faultCountIndex=faultIdx390  % cycle over each fault
    % linear interpolation to find pressure to failture at desired
    % probability
   PressureToSlipAtQueriedProbability(faultCountIndex,1) = interp1(cdfLineDataThisAreaYData(2:end,faultCountIndex),cdfLineDataThisAreaXData(2:end,faultCountIndex),probabilityToQueryAt,'linear' ); %
   
end

% draw red x
scatter(PressureToSlipAtQueriedProbability,ones(size(PressureToSlipAtQueriedProbability)).*probabilityToQueryAt,'rx')




