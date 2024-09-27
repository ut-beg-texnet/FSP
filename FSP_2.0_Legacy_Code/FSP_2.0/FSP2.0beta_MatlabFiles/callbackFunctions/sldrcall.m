


% year slider callback
function sldrcall(src,~,hDV)

val = get(src,'value') ; val = floor(val) ;

if hDV.data.reservoir.importHydrology==1 % what if hydrologic model is imported and not all years are contiguous?
    if ~any(hDV.data.reservoir.yearsRepresentedHydroImport==val) % see if a model was imported for that year
        % then find next lowest value
        sotredYears=sort([hDV.data.reservoir.yearsRepresentedHydroImport;val]);% sort the desired year in years available
        idx=find(sotredYears==val); % find where this year is in loaded years
        if idx~=length(sotredYears)-1 % if not last entry
            val=sotredYears(idx+1); % go one model later
        else
            val=sotredYears(idx-1);% go one model earlier
        end
    end
end
mival = get(src,'min') ; mxval = get(src,'max') ;
if val<mival
    val=mival;
elseif val>mxval
    val=mxval;
end
set(src,'value',val) ;
%         set(yrtxt,'string',num2str(val)) ;

%update all other sliders and text
set(hDV.hdsldr(:),'value',val) ;
set(hDV.hdsldr_txt(:),'string',num2str(val)) ;

if hDV.plotdata.printFunctionName
    disp(['running slider callback for ',num2str(val),' about to call calcengine and refreshplotdata'])
end
calcengine(hDV,hDV.currtab.name); %rerun calculations
refreshplotdata(hDV,hDV.currtab.name); %plot data for selected panel refreshed
end


