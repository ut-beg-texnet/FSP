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
% inj inputted variables

% modified by Rall to save inputted variables too

function [out, inj]= monte_carlo(f,in0,inSig,hDV,nruns)

% disp('using monte_carloSaveInsDeleteme, not monte carlo ')
out = zeros(nruns,1);% preallocate
inj=cell(nruns,length(in0));% preallocate
for jj=1:1:nruns
 
%     inj=in0; % reassign
    for k=1:1:length(in0)

        %inj{k} = inj{k}+inSig{k}.*randn(1,length(inSig{k})); % use this for normal distribution
        inj{jj,k} = in0{k}+inSig{k}.*(2*rand(1,length(inSig{k}))-1); % use this for uniform distribution
    end
    out(jj) = f(inj(jj,:),hDV) ;
end

end