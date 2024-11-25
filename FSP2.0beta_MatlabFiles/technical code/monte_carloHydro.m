% D. Pais
% ExxonMobil Upstream Reseach Company 2016
% Basic Monte Carlo Engine

% Inputs: f: function handle of map to run Monte Carlo On
%         in0: cell array of nominal inputs to function, i.e. f(in0)
%         produces the function output
%         inSig: sigma (gaussian) cell array for the sigma associated with
%         each element of in0 ; or delta +/- in uniform distribution
%         nruns: how many Monte Carlo Runs you want to do

% Out: desired outputs

% Modified by Rall to ignore but pass through cell inputs, which can't be
% added to

function [out,inj] = monte_carloHydro(f,in0,inSig,nruns,nAnswers)
if nargin==4
    nAnswers=1;
end
        
mainFigure=gcf;
%wait bar
h = waitbar(0.01,'Calculating ...','Name','Monte-Carlo Engine');
centerFigure(mainFigure,h)

out = zeros(nruns,nAnswers);
inj=cell(nruns,length(in0));% preallocate
for j=1:1:nruns
    
    %     inj=in0; % reassign
    for k=1:1:length(in0)
        if ~iscell(in0{k}) % ignore cell data of well data, ie. allWellsDatenumBarrelsPerDay
            %inj{k} = inj{k}+inSig{k}.*randn(1,length(inSig{k})); % use this for normal distribution
            inj{j,k} = in0{k}+inSig{k}.*(2*rand(1,length(inSig{k}))-1); % use this for uniform distribution
        else % cell within a cell - don't vary injection data
            inj{j,k} = in0{k};
        end
    end
    out(j,:) = f(inj(j,:))' ;
    
    if ~ishandle(h)
        msgWindow1=msgbox('Monte Carlo Premature Shutdown', 'Monte Carlo Premature Shutdown warning','warn');
        centerFigure(hDV.hfig,msgWindow1);
        break
    end
      if mod(j,5)==1 % don't update every run of loop, faster
        waitbar(j / nruns,h,['Calculating Hydrology ',num2str(100*j / nruns,'%3.0f'), '% ...']) ;
     end
end

%remove wait bar
if ishandle(h)
    close(h) ;
end


end