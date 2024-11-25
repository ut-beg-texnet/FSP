% by Rall Walsh
% stanford Summer 2016
% read in text data of disposal well location and save in data format that
% can be plotted: MATLAB datenumber and barrelsPerDay
% currently assumes the rate in the last monthly entry continues till the end of the
% model
% this function could be improved by checking for repeat monthly data
% entries.
% note that datapoint plots at last day with that injection rate, not first
% day with that injection rate. 

function SpreadsheetStrings2WellData(stringsWellDataAdvanced,hDV)

%%%% only for testing
if nargin ==0 % for coding
    load('test code/testData.mat')
    hDV.data.nwells_max=25;
    hDV.data.realWellData.columnIsNumber=columnIsNumber;
    hDV.data.realWellData.stringsWellDataAdvanced=stringsWellDataAdvanced;
    hDV.data.realWellData.extrapolateInjectionCheck=1;
    hDV.hfig=figure;
    %       nargin=0;
end
%%%%


numberColumns=hDV.data.realWellData.stringsWellDataAdvanced(:,hDV.data.realWellData.columnIsNumber); % find numerical data
numbers=str2double(numberColumns); % convert numerical data from strings to numbers

uniqueWellColumnIndex=1; % find unique wells based on unique entries of column 1: API/name
[C,IA,IC]=unique(hDV.data.realWellData.stringsWellDataAdvanced(:,uniqueWellColumnIndex));
hDV.data.realWellData.wellNames=hDV.data.realWellData.stringsWellDataAdvanced(IA,1);% well name stored once per well
% hDV.data.realWellData.wellTypes=hDV.data.realWellData.stringsWellDataAdvanced(IA,9);% well type stored once per well
% hDV.data.realWellData.XEasting=str2double(hDV.data.realWellData.stringsWellDataAdvanced(IA,2));% well X stored once per well
% hDV.data.realWellData.YNorthing=str2double(hDV.data.realWellData.stringsWellDataAdvanced(IA,3));% well Y stored once per well
hDV.data.nwells=size(C,1); % count number of wells

if  hDV.data.nwells > hDV.data.nwells_max % if more wells than max number
    %     waitfor(msgbox(cat(2,'Currently FSP can handle ',num2str(hDV.data.nwells_max),' wells, but based on unique API number/name column, we see ',num2str(hDV.data.nwells),...
    %         ' wells. Make sure you aren''t entering the same well with different entries here. Just taking the first ',num2str(hDV.data.nwells_max),' wells and ignoring the rest.'),'csv load','help'));
    errorWindow1= errordlg(cat(2,'Currently FSP can handle ',num2str(hDV.data.nwells_max),' wells, but based on unique API number/name column, we see ',num2str(hDV.data.nwells),...
        ' wells. Make sure you aren''t entering the same well with different entries here. Just taking the first ',num2str(hDV.data.nwells_max),' wells and ignoring the rest.'));
    centerFigure(hDV.hfig,errorWindow1);
    hDV.data.nwells=hDV.data.nwells_max; % then only load in nwellsMax
    % and throw out other data
%     hDV.data.realWellData.XEasting=hDV.data.realWellData.XEasting(1:hDV.data.nwells_max);
%     hDV.data.realWellData.YNorthing=hDV.data.realWellData.YNorthing(1:hDV.data.nwells_max);
    hDV.data.realWellData.wellNames=hDV.data.realWellData.wellNames(1:hDV.data.nwells_max);
else
    msgWindow1= msgbox(cat(2,'Based on unique API number/name, we see ',num2str(hDV.data.nwells),...
        ' wells in the spreadsheet. Check that column if you think this is incorrect.'),'csv load','help');
    centerFigure(hDV.hfig,msgWindow1);
end

for CycleWellsK=1:hDV.data.nwells % cycle over wells
     
    % warn if all lat/longs are not the same for a given unique well name
     allEastAndNorthThisWell=str2double(hDV.data.realWellData.stringsWellDataAdvanced(IC==CycleWellsK,2:3));
     UniqueEastNorthThisWell=unique(allEastAndNorthThisWell,'rows');
     hDV.data.realWellData.YNorthing(CycleWellsK,1)=allEastAndNorthThisWell(1,2);
     hDV.data.realWellData.XEasting(CycleWellsK,1)=allEastAndNorthThisWell(1,1);

     % if not all X oordinates are the same for this well
    if     size(UniqueEastNorthThisWell,1)>1 % if more than 1 set of coordinates 
        coordChosenRepeated=[ones(size(UniqueEastNorthThisWell,1),1).*hDV.data.realWellData.XEasting(CycleWellsK,1),ones(size(UniqueEastNorthThisWell,1),1).*hDV.data.realWellData.YNorthing(CycleWellsK,1)];
        difflocs=UniqueEastNorthThisWell(:,1)+1i*UniqueEastNorthThisWell(:,2) - (coordChosenRepeated(:,1)+1i*coordChosenRepeated(:,2)) ; distancesBetweenCoordiantes=abs(difflocs) ;    %grid distances to well
        msgWindow1=msgbox({cat(2,'Warning, not all East/North coordinates for well: ',hDV.data.realWellData.wellNames{CycleWellsK},' are the same. Coordinates were:');...
            'East,            North,            meters from chosen coordinate:';...            
            num2str([UniqueEastNorthThisWell,distancesBetweenCoordiantes.*1000]);...
            cat(2,'but the coordinates ',num2str(hDV.data.realWellData.XEasting(CycleWellsK)),', ',num2str(hDV.data.realWellData.YNorthing(CycleWellsK)),' was arbitrarily chosen')},...
            'Well coordinate inconsistency warning','warn');
        centerFigure(hDV.hfig,msgWindow1);
    end
    
    % these will be hDV.data.realWellData. ...
    %     hDV.data.realWellData.APIs(CycleWellsK,1)=C(CycleWellsK,1);
    %     hDV.data.realWellData.XEasting(CycleWellsK,1)=C(CycleWellsK,2);
    %     hDV.data.realWellData.YNorthing(CycleWellsK,1)=C(CycleWellsK,3);
    
    
    % should make it identify repeated monthly data here, and throw warning/throw
    % out data.
    YearMonthVolumeThisWell=numbers(IC==CycleWellsK,3:5); % find data for this well
    YearMonthVolumeThisWell=sortrows(YearMonthVolumeThisWell,[1,2])   ;% sort it sequentially based on year and month
    
    hDV.data.realWellData.YearMonthVolume{CycleWellsK,1}=YearMonthVolumeThisWell; % save this well numerical data
    
    datenumBarrelsPerDayThisWell=[];  % empty
    
    % cycle over each month in this well's data
    for CycleMonthsK=1:size(YearMonthVolumeThisWell,1)
        
        % find matlab datenum of year and month
        thisYear=YearMonthVolumeThisWell(CycleMonthsK,1);
        thisMonth=YearMonthVolumeThisWell(CycleMonthsK,2);
        barrelsPerDayOverThisMonth=YearMonthVolumeThisWell(CycleMonthsK,3)./eomday(thisYear,thisMonth); % divide monthly volume by number of days in Month
        %         startDateThisMonthdatenum=datenum(thisYear,thisMonth,1,0,0,0); % beginning of first day of month
        endDateThisMonthdatenum=datenum(thisYear,thisMonth,eomday(thisYear,thisMonth),23,59,59); % last second of this month
        
        % find previous month
        previousMonth=thisMonth-1;
        if previousMonth==0 % if this Month is January
            previousMonth=12; % previous month is December
            previousYear=thisYear-1; % year of previous month
        else
            previousYear=thisYear;  % year of previous month
        end
        
        % if there is a data entry for previous year and month:
        if any(YearMonthVolumeThisWell(:,1)==previousYear & YearMonthVolumeThisWell(:,2)==previousMonth)
            %  then do nothing because that data entry exists
            
            % if injection rate didn't change from last month to this month:
            if YearMonthVolumeThisWell(YearMonthVolumeThisWell(:,1)==previousYear & YearMonthVolumeThisWell(:,2)==previousMonth,3)==...
                    YearMonthVolumeThisWell(YearMonthVolumeThisWell(:,1)==thisYear & YearMonthVolumeThisWell(:,2)==thisMonth,3)
                % then last month's injection rate equals this months injection rate, so
                % move that time data point from the end of last month to the
                % end of this month.% this makes the datapoints not duplicate
                % when wells injct at a constant rate.
                % this doesn't account for the number of days in each month, unlike the
                % other parts of the code.
                datenumBarrelsPerDayThisWell(end,1)=endDateThisMonthdatenum; % move last months datenumber to this months datenumber, rate doesn't change
                continue % go to next month in for loop, just moved time point from previous month.
            end
        else % there is no entry for previous month, so make last second of previous month data zero injection rate
            endDatePreviousMonthdatenum=datenum(previousYear,previousMonth,eomday(previousYear,previousMonth),23,59,59); % last second of previous month
            datenumBarrelsPerDayThisWell(end+1,[1:2])=[endDatePreviousMonthdatenum,0]; % make entry of 0 barrels previous month
        end
        
        % this is the line that matters the most:
        datenumBarrelsPerDayThisWell(end+1,[1:2])=[endDateThisMonthdatenum,barrelsPerDayOverThisMonth]; % make entry at end of month for this month
        
    end % done cycling over months of this well
    
    
    
    %     assumeInjectionExtendsToEndOfModel=true; % could make this a checkbox
    if hDV.data.realWellData.extrapolateInjectionCheck % checkbox for extrapolation value
        if nargin ==0 % for coding
            YearEndOfModel=2045;
        else
            YearEndOfModel=get(hDV.hdsldr(1),'max');% end of scaleruler.
        end
        %
        datenumBarrelsPerDayThisWell(end,1)=[datenum(YearEndOfModel,1,1,0,1,0)];% move last point to last year in model
    end
    
    hDV.data.realWellData.datenumBarrelsPerDay{CycleWellsK,1}=datenumBarrelsPerDayThisWell; % save datenumber and bblsPerDay
    
end % done cycling over all wells

% Clear any data of previous wells beyond hDV.data.nwells:
% hDV.data.realWellData.APIs=hDV.data.realWellData.APIs(1:hDV.data.nwells,1);
hDV.data.realWellData.XEasting=hDV.data.realWellData.XEasting(1:hDV.data.nwells,1);
hDV.data.realWellData.YNorthing=hDV.data.realWellData.YNorthing(1:hDV.data.nwells,1);
hDV.data.realWellData.wellNames=hDV.data.realWellData.wellNames(1:hDV.data.nwells,1);
% hDV.data.realWellData.wellTypes=hDV.data.realWellData.wellTypes(1:hDV.data.nwells,1);
hDV.data.realWellData.YearMonthVolume=hDV.data.realWellData.YearMonthVolume(1:hDV.data.nwells,1);
hDV.data.realWellData.datenumBarrelsPerDay=hDV.data.realWellData.datenumBarrelsPerDay(1:hDV.data.nwells,1);

if nargin==0
    % check by plotting
    
    hold on
    for dd=1:length(hDV.data.realWellData.datenumBarrelsPerDay)
        datenumBarrelsPerDayThisWell=hDV.data.realWellData.datenumBarrelsPerDay{dd,1};
        [tt,QQ]=stairs(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2));
        plot(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2),'rx') % data points as red xs
        
        % do this adjustment because value represents rate up to that point, not
        % rate after that point as stair assumes.
        % see plottingInjectionStairData.m
        QQ(1:2)=[];
        QQ(end+1:end+2,1)=QQ(end);
        % these points bring the line down to the X axis at the beginning and end
        QQ=[0;QQ;0];
        tt=[tt(1);tt;tt(end)];
        plot(tt,QQ,'b') % injection pattern in blue
    end
    
    datetick('x')
    
    
end % end SpreadsheetStrings2WellData function











