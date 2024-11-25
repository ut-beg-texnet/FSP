% check probabilistic data distribution bound entries
% second entry was noshow, but don't care about here
%
% Modified to check only if probabilistic is to be used, especially not in
% case hydrology is imported from reservoir model ...
% S P Lele, Oct 2018
%
function [flag] = checkProbabilisticEntriesHyd(hDV,noshow,vec,txt)

if hDV.plotdata.printFunctionName
    disp(['running checkProbabilisticEntriesHyd '])
end

flag=true; % flag=1 if no issues
% by Rall
% still need to check relative stress magnitudes etc, and if stress state
% is possible.
deterministicVals=[hDV.data.reservoir.vals(1);hDV.data.reservoir.vals(2);hDV.data.reservoir.vals(3);hDV.data.adv.vals(5);hDV.data.adv.vals(7);hDV.data.adv.vals(8);hDV.data.adv.vals(9);0];

switch nargin % number of function entries
    case 3
        txt=hDV.data.probHydrology.distirbutionTxt;
    case 2
        txt=hDV.data.probHydrology.distirbutionTxt;
        if isfield(hDV.data.probHydrology, 'sigvals')
            vec = hDV.data.probHydrology.sigvals;
        else
            vec = zeros(1,8);
        end
    case 1
        txt=hDV.data.probHydrology.distirbutionTxt;
        if isfield(hDV.data.probHydrology, 'sigvals')
            vec = hDV.data.probHydrology.sigvals;
        else
            vec = zeros(1,8);
        end
        noshow=0;
end

% Run checks only if probabilistic to be used
%
if hDV.data.probHydrology.probabilistic1VsDeterministic2==1
    for kj=1:1:8;
        if vec(kj)<0
            edlgbox = errordlg(cat(2,'check ',txt{kj},' you entered a negative number: ',num2str(vec(kj)),...
                ' FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution'));
            centerFigure(hDV.hfig,edlgbox);
            flag=false;
        elseif vec(kj)>=deterministicVals(kj) && kj~=8 % ignore bootstraps for this criteria
            edlgbox = errordlg(cat(2,'check ',txt{kj},' you entered an uncertainty: ',num2str(vec(kj)),...
                ' greater than or equal to the deterministic value. FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution',...
                ' and you can''t have a negative value for this parameter'));
            centerFigure(hDV.hfig,edlgbox);
            flag=false;
        end
        if kj==2 && deterministicVals(kj)+vec(kj)>=100; % porosity less than 1
            edlgbox = errordlg(cat(2,'check ',txt{kj},' you entered an uncertainty: ',num2str(vec(kj)),...
                ' That allows poroity to be 100% or more. FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution',...
                ' and you can''t have a porosity >=100%'));
            centerFigure(hDV.hfig,edlgbox);
            flag=false;
        end

        if kj==8 && vec(kj)<1  || kj==8 && vec(kj)~=round(vec(kj))  ; % bootstrap number
            edlgbox = errordlg(cat(2,'check number of hydrologic iterations. You entered a number of iterations: ',num2str(vec(kj)),...
                ' That is not a positive integer. This is the number of times the hydrologic calculation is run.',...
                ' Higher numbers take longer to calculate, but give smoother probability of exceedance curves. '));
            centerFigure(hDV.hfig,edlgbox);
            flag=false;
        end
    end
end
end
