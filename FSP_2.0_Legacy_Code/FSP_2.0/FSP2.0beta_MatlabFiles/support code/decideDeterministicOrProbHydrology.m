% by Rall
% Stanford 2016
% decide between deterministic and probabilistic hydrology
% if value is already set, keep it.
% can't do probabilistic hydrology with an imported modflow model
% Should do probabilistic hydrology with constant rate injection
% optional with monthly injection data, depending on the amount, as lots of
% months runs slowly
%
% Moved checks whether realWellData used etc. to block corresponding to
% importHydrology==0; if hydrology imported then well data not used
% And ef the variable probabilistic1VsDeterministic2 already exists, still
% add check for importHydrology instead of returning value as is
% --SPL, Sept 2018
%
function probabilistic1VsDeterministic2= decideDeterministicOrProbHydrology(hDV)

if ~isfield(hDV.data.probHydrology,'probabilistic1VsDeterministic2') % if variable doesn't exist
        
    if hDV.data.reservoir.importHydrology   % imported modflow model? 1=yes
        % can't do probabilistic hydrology on a modflow model
        probabilistic1VsDeterministic2=2;
    else
        if isfield(hDV.data.realWellData, 'use') && hDV.data.realWellData.use  %if variable rate injection data loaded
            % this is the marginal/optional/decision case
            % can be slow, but also can be desireable

            % so if lots of monthly data, default to  deterministic:
            % arbitrarily chose 400 months of data s cutoff
            % decide how many bootstraps based on data size
            if size(hDV.data.realWellData.stringsWellDataAdvanced,1)>400
                probabilistic1VsDeterministic2=2;
            elseif size(hDV.data.realWellData.stringsWellDataAdvanced,1)>200
                probabilistic1VsDeterministic2=1;
                hDV.data.probHydrology.sigvals(8)=400; % set 400 bootstraps
            else % few well-months, so can do more bootstraps
                % if less than that amount of months of data, default to probabilistic
                probabilistic1VsDeterministic2=1;
                hDV.data.probHydrology.sigvals(8)=1000; % set 1000 bootstraps
            end
        else % make constant rate injection into format of variable rate data
            probabilistic1VsDeterministic2=1;
        end
    end
    % assign whatever was chosen
    hDV.data.probHydrology.probabilistic1VsDeterministic2=probabilistic1VsDeterministic2;
    
else % if value already exists
    
    if hDV.data.reservoir.importHydrology   % imported modflow model? 1=yes
        % can't do probabilistic hydrology on a modflow model
        probabilistic1VsDeterministic2=2;
    else
        probabilistic1VsDeterministic2=hDV.data.probHydrology.probabilistic1VsDeterministic2;
    end
end