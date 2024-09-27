% by Rall Walsh
% stanford Summer 2016
% read in text data of disposal well location and save in data format that
% can be plotted: decyear and barrelsPerDay
% assume the rate in the last monthly entry continues til the end of the
% model
% also save in version that can calculate pressure from variuable rate
% well: datenum and bbl/day

function SpreadsheetStrings2WellDataOld2(stringsWellDataAdvanced,hDV)

%%%% only for testing
if nargin ==0 % for coding
load('testData.mat')
hDV.data.nwells_max=100;
hDV.data.realWellData.columnIsNumber=columnIsNumber;
hDV.data.realWellData.stringsWellDataAdvanced=stringsWellDataAdvanced;
end
%%%%


numberColumns=hDV.data.realWellData.stringsWellDataAdvanced(:,hDV.data.realWellData.columnIsNumber); % find numerical data
numbers=str2double(numberColumns); % convert numerical data to numbers

uniqueWellColumnIndexes=1; % find unique wells based on unique entries of rows 1:3, API, east and north
[C,IA,IC]=unique(numbers(:,uniqueWellColumnIndexes),'rows');
hDV.data.realWellData.wellNames=hDV.data.realWellData.stringsWellDataAdvanced(IA,8);% well name stored once per well
hDV.data.realWellData.wellTypes=hDV.data.realWellData.stringsWellDataAdvanced(IA,9);% well type stored once per well
hDV.data.realWellData.XEasting=str2double(hDV.data.realWellData.stringsWellDataAdvanced(IA,2));% well X stored once per well
hDV.data.realWellData.YNorthing=str2double(hDV.data.realWellData.stringsWellDataAdvanced(IA,3));% well Y stored once per well

hDV.data.nwells=size(C,1); % count number of wells

if  hDV.data.nwells > hDV.data.nwells_max % if more wells than max number
    waitfor(msgbox(cat(2,'Currently FSP can handle ',num2str(hDV.data.nwells_max),' wells, but based on unique Easting, Northing, and API numbers, we see ',num2str(hDV.data.nwells),...
        ' wells. Make sure you aren''t entering the same well with different API or coordinates. Just taking the first ',num2str(hDV.data.nwells_max),' wells and ignoring the rest.'),'csv load','help'));
        hDV.data.nwells=hDV.data.nwells_max; % then only load in nwellsMax
else
    waitfor(msgbox(cat(2,'Based on unique Easting, Northing, and API numbers, we see ',num2str(hDV.data.nwells),...
        ' wells in the spreadsheet. Check API numbers if you think this is incorrect.'),'csv load','help'));
end


for CycleWellsK=1:hDV.data.nwells % cycle over wells
    
    % these will be hSV.data.realWellData. ...
    hDV.data.realWellData.APIs(CycleWellsK,1)=C(CycleWellsK,1);
%     hDV.data.realWellData.XEasting(CycleWellsK,1)=C(CycleWellsK,2);
%     hDV.data.realWellData.YNorthing(CycleWellsK,1)=C(CycleWellsK,3);
    

% should make it identify repeated monthly data here, and throw warning/throw
% out data. 
    YearMonthVolumePressureThisWell=numbers(IC==CycleWellsK,4:7); % find data for this well
    YearMonthVolumePressureThisWell=sortrows(YearMonthVolumePressureThisWell,[1,2])   ;% sort it sequentially
    
    hDV.data.realWellData.YearMonthVolumePressure{CycleWellsK,1}=YearMonthVolumePressureThisWell; % save this well numerical data
    
    % clear
    decyearBarrelsPerDayThisWell=[];
    datenumBarrelsPerDayThisWell=[];  % this is for a stair plot 

    % cycle over each month in this well's data
    for CycleMonthsK=1:size(YearMonthVolumePressureThisWell,1)
        
%                 volumeOfPreviousMonth=[]; % clear out 
                
                
        % find matlab decyear of year and month
        thisYear=YearMonthVolumePressureThisWell(CycleMonthsK,1);
        thisMonth=YearMonthVolumePressureThisWell(CycleMonthsK,2);
        startDateThisMonthdecyear=decyear(thisYear,thisMonth,1,0,0,0); % beginning of first day of month
        endDateThisMonthdecyear=decyear(thisYear,thisMonth,eomday(thisYear,thisMonth),23,59,59); % last second of month
        barrelsPerDayOverThisMonth=YearMonthVolumePressureThisWell(CycleMonthsK,3)./eomday(thisYear,thisMonth); % divide monthly volume by number of days in Month
        startDateThisMonthdatenum=datenum(thisYear,thisMonth,1,0,0,0); % beginning of first day of month
%         endDateThisMonthdatenum=decyear(thisYear,thisMonth,eomday(thisYear,thisMonth),23,59,59); % last second of month
        
        
        % find previous month
        previousMonth=thisMonth-1;
        if previousMonth==0 % if this Month is January
            previousMonth=12; % previous month is December
            previousYear=thisYear-1; % year of previous month
        else
            previousYear=thisYear;  % year of previous month
        end
        
        % find month after
        nextMonth=thisMonth+1;
        if nextMonth==13 % if this Month is December
            nextMonth=1; % next month is January
            nextYear=thisYear+1; % year of previous month
        else
            nextYear=thisYear;  % year of next month
        end
        
        % find 2 months after this month
        twomonthsAfter=thisMonth+2;
        if twomonthsAfter==13 % if this Month is November
            twomonthsAfter=1; % in 2 months, it's January
            yearInTwoMonths=thisYear+1; % year in 2 months
        elseif twomonthsAfter==14 % if this Month is December
            twomonthsAfter=2; %  month is feb
            yearInTwoMonths=thisYear+1; % year in two months    
        else
            yearInTwoMonths=thisYear;  % year in two months 
        end
        
        
        
        % if there is a data entry for previous year and month
        if any(YearMonthVolumePressureThisWell(:,1)==previousYear & YearMonthVolumePressureThisWell(:,2)==previousMonth) 
            % do nothing because that data entry exists
%             volumeOfPreviousMonth=YearMonthVolumePressureThisWell(YearMonthVolumePressureThisWell(:,1)==previousYear & YearMonthVolumePressureThisWell(:,2)==previousMonth,3);
        else % no entry for previous month, so make it zero
            endDatePreviousMonthdecyear=decyear(previousYear,previousMonth,eomday(previousYear,previousMonth),23,59,59); % last second of previous month
            decyearBarrelsPerDayThisWell(end+1,[1:2])=[endDatePreviousMonthdecyear,0]; % make entry of 0 barrels previous month
            startDatePreviousMonthdatenum=datenum(previousYear,previousMonth,1,0,0,0); % first second of previous month         
           datenumBarrelsPerDayThisWell(end+1,[1:2])=[startDatePreviousMonthdatenum,0]; % make entry of 0 barrels previous month 
        end
        
                decyearBarrelsPerDayThisWell(end+1,[1:2])=[startDateThisMonthdecyear,barrelsPerDayOverThisMonth]; % make entry at beginning of month  (green dot)
 

        
        decyearBarrelsPerDayThisWell(end+1,[1:2])=[endDateThisMonthdecyear,barrelsPerDayOverThisMonth]; % make entry at end of month  (green dot)
 
        datenumBarrelsPerDayThisWell(end+1,[1:2])=[startDateThisMonthdatenum,barrelsPerDayOverThisMonth]; % make entry at beginning of month  (black dot)
 
     
        if any(YearMonthVolumePressureThisWell(:,1)==nextYear & YearMonthVolumePressureThisWell(:,2)==nextMonth) % if there is a data entry for next month
            % do nothing because that data entry exists
        else % no entry for next month starting with 0
            startDateNextMonthdecyear=decyear(nextYear,nextMonth,1,0,0,0); % beginning of first day of next month
            decyearBarrelsPerDayThisWell(end+1,[1:2])=[startDateNextMonthdecyear,0]; % make entry of 0 barrels next month (green dot)
            
            if ~any(YearMonthVolumePressureThisWell(:,1)==yearInTwoMonths & YearMonthVolumePressureThisWell(:,2)==twomonthsAfter) % if there is no data entry for two months after
            startDateNextMonthdatenum=datenum(nextYear,nextMonth,1,0,0,0); % beginning of first day of next month       
            datenumBarrelsPerDayThisWell(end+1,[1:2])=[startDateNextMonthdatenum,0]; % make entry at beginning of month  (black x)
            end
            
        end
        

    end % done cycling over months of this well
    
    
    
    assumeInjectionExtendsToEndOfModel=1;
    
    if assumeInjectionExtendsToEndOfModel
        maxyear=2045; % need a way to find this from max(slider years) 
        decyearBarrelsPerDayThisWell(end,2)=[decyearBarrelsPerDayThisWell(end-1,2)]; % last second of last year in model
        decyearBarrelsPerDayThisWell(end+1,1:2)=[decyear(maxyear-1,12,31,23,59,59),decyearBarrelsPerDayThisWell(end,2)]; % last second of last year in model
        
        datenumBarrelsPerDayThisWell(end,2)=[datenumBarrelsPerDayThisWell(end-1,2)]; 
        datenumBarrelsPerDayThisWell(end+1,1:2)=[datenum(maxyear-1,12,31,23,59,59),datenumBarrelsPerDayThisWell(end,2)]; % last second of last year in model
       
        
    end
    
    hDV.data.realWellData.decyearBarrelsPerDay{CycleWellsK,1}=decyearBarrelsPerDayThisWell; % save decyear and bblsPerDay that will be plotted
    hDV.data.realWellData.datenumBarrelsPerDay{CycleWellsK,1}=datenumBarrelsPerDayThisWell; % save datenumber and bblsPerDay , you need to do stairs(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2)) to plot this
    
    
end % done cycling over all wells



% Clear any data of previous wells beyond hDV.data.nwells:
hDV.data.realWellData.APIs=hDV.data.realWellData.APIs(1:hDV.data.nwells,1);
hDV.data.realWellData.XEasting=hDV.data.realWellData.XEasting(1:hDV.data.nwells,1);
hDV.data.realWellData.YNorthing=hDV.data.realWellData.YNorthing(1:hDV.data.nwells,1);
hDV.data.realWellData.wellNames=hDV.data.realWellData.wellNames(1:hDV.data.nwells,1);
hDV.data.realWellData.wellTypes=hDV.data.realWellData.wellTypes(1:hDV.data.nwells,1);
hDV.data.realWellData.YearMonthVolumePressure=hDV.data.realWellData.YearMonthVolumePressure(1:hDV.data.nwells,1);
hDV.data.realWellData.decyearBarrelsPerDay=hDV.data.realWellData.decyearBarrelsPerDay(1:hDV.data.nwells,1);
hDV.data.realWellData.datenumBarrelsPerDay=hDV.data.realWellData.datenumBarrelsPerDay(1:hDV.data.nwells,1);


% d=hDV.data.realWellData.decyearBarrelsPerDay




if nargin==0
    % check by plotting 
   % plot 
   figure
   hold on
   for dd=1:length(hDV.data.realWellData.decyearBarrelsPerDay)
   plot(hDV.data.realWellData.decyearBarrelsPerDay{dd}(:,1),hDV.data.realWellData.decyearBarrelsPerDay{dd}(:,2),'linewidth',3)
   plot(hDV.data.realWellData.decyearBarrelsPerDay{dd}(:,1),hDV.data.realWellData.decyearBarrelsPerDay{dd}(:,2),'go','markersize',6,'linewidth',3)   
   
   end
%    figure
   for dd=1:length(hDV.data.realWellData.decyearBarrelsPerDay)
   stairs(decyear(datevec(hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,1))),hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,2),'r:','linewidth',2)
   plot(decyear(datevec(hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,1))),hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,2),'kx','linewidth',3,'markersize',3)
   [x6,y6]=  stairs(decyear(datevec(hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,1))),hDV.data.realWellData.datenumBarrelsPerDay{dd}(:,2))
   plot(x6,y6,'ro','markersize',12);
   end


end % end SpreadsheetStrings2WellData function











