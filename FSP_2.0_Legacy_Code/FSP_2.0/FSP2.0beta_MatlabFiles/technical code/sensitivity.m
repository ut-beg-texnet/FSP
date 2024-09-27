% D. Pais
% ExxonMobil Upstream Reseach Company 2016
% Basic Sensitivity Analysis

% Inputs: f: function handle of map to run sensitivity On
%         in0: cell array of nominal inputs to function, i.e. f(in0)
%         produces the function output
%         inSig: sigma (gaussian) cell array for the sigma associated with
%         each element of in0
%         ptile: percentile level

% Outut: outup, outdown: outputs for increasing/decreasing the input parameters
%        nom: nominal output (f(in0))

function [outup , outdown , nom] = sensitivity(f,in0,inSig,hDV,ptile)

if nargin==4
    ptile=0.001; 
end

%placeholders
outup=in0;
outdown=in0;

%nominal value
nom=f(in0,hDV);

for k=1:1:length(in0)
    valup=nom+0*in0{k} ; valdown=nom+0*in0{k}; %paceholders for outputs
    temp=in0{k}; tempsig=inSig{k};
    for i=1:1:length(in0{k})
        if tempsig(i)~=0
            inup=in0; indown=in0;
            
            %[upval,downval]=getptilevals(temp(i),tempsig(i),ptile); %use this for normal distribution
            
            upval = temp(i)+tempsig(i) ; downval = temp(i)-tempsig(i) ; %use this for uniform distribution full range
            
            tup=temp;tup(i)=upval;
            tdown=temp; tdown(i)=downval;
            inup{k}=tup;
            indown{k}=tdown;
            valup(i)=f(inup,hDV);
            valdown(i)=f(indown,hDV);
        end
    end
    outup{k} = valup;
    outdown{k} = valdown;
end

    % Function to be used if normal distributions for input variables are
    % used
    %
    function [up down] = getptilevals(mu,sig,ptile)
        up = norminv(1-ptile,mu,sig);
        down = norminv(ptile,mu,sig);
    end



end
