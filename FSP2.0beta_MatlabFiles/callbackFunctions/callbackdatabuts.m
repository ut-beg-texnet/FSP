%call back function on clicking data buttons
%this must tie to butcall in dataentry.m
function callbackdatabuts(src,event,i,hDV)

switch i
    case 1
        dataentryStress(hDV,get(src,'Label'),hDV.data.stress.txt,hDV.data.stress.vals,i);
    case 2
        dataentryHydrology(hDV,get(src,'Label'),hDV.data.reservoir.txt,hDV.data.reservoir.vals,i,[4 5]);
    case 3
            if   hDV.data.reservoir.importHydrology==0 % if you load a hydrologic pressure model, you don't need well data
                dataentryWells(hDV,get(src,'Label'),hDV.data.inject.txt,hDV.data.inject.vals,i);
            else
                msgWindow1=msgbox({cat(2,'Under "Hydrology Data" you selected loading an external hydrologic pressure model, so you can''t enter wells')}, 'Loading Hydrologic Model','warn');
                centerFigure(hDV.hfig,msgWindow1);
            end
    case 4
        dataentryFaults(hDV,get(src,'Label'),hDV.data.fault.txt,hDV.data.fault.vals,i);
        
    case 5
        dataentryAdvanced(hDV,get(src,'Label'),hDV.data.adv.txt,hDV.data.adv.vals,i,[6]);
        
end


end

