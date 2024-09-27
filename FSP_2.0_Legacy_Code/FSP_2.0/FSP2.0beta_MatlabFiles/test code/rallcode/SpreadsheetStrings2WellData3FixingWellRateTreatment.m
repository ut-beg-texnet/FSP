% by Rall Walsh
% stanford Summer 2016
% read in text data of disposal well location and save in data format that
% can be plotted: decyear and barrelsPerDay
% assume the rate in the last monthly entry continues til the end of the
% model
% also save in version that can calculate pressure from variuable rate
% well: datenum and bbl/day

function hDV=SpreadsheetStrings2WellData3FixingWellRateTreatment(stringsWellDataAdvanced,hDV)

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
    waitfor(msgbox(cat(2,'Currently FSP can handle ',num2str(hDV.data.nwells_max),' wells, but based on unique API numbers, we see ',num2str(hDV.data.nwells),...
        ' wells. Make sure you aren''t entering the same well with different API. Just taking the first ',num2str(hDV.data.nwells_max),' wells and ignoring the rest.'),'csv load','help'));
    hDV.data.nwells=hDV.data.nwells_max; % then only load in nwellsMax
else
    waitfor(msgbox(cat(2,'Based on unique API numbers, we see ',num2str(hDV.data.nwells),...
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
    YearMonthVolumePressureThisWell=sortrows(YearMonthVolumePressureThisWell,[1,2])   ;% sort it sequentially based on year and month
    
    hDV.data.realWellData.YearMonthVolumePressure{CycleWellsK,1}=YearMonthVolumePressureThisWell; % save this well numerical data
    
    % clear
    %     decyearBarrelsPerDayThisWell=[];
    datenumBarrelsPerDayThisWell=[];  % this is for a stair plot
    
    % cycle over each month in this well's data
    for CycleMonthsK=1:size(YearMonthVolumePressureThisWell,1)
        
        %                 volumeOfPreviousMonth=[]; % clear out
        
        
        % find matlab decyear of year and month
        thisYear=YearMonthVolumePressureThisWell(CycleMonthsK,1);
        thisMonth=YearMonthVolumePressureThisWell(CycleMonthsK,2);
        %         startDateThisMonthdecyear=decyear(thisYear,thisMonth,1,0,0,0); % beginning of first day of month
        %         endDateThisMonthdecyear=decyear(thisYear,thisMonth,eomday(thisYear,thisMonth),23,59,59); % last second of month
        barrelsPerDayOverThisMonth=YearMonthVolumePressureThisWell(CycleMonthsK,3)./eomday(thisYear,thisMonth); % divide monthly volume by number of days in Month
        startDateThisMonthdatenum=datenum(thisYear,thisMonth,1,0,0,0); % beginning of first day of month
        endDateThisMonthdatenum=datenum(thisYear,thisMonth,eomday(thisYear,thisMonth),23,59,59); % last second of this month
        
        % find previous month
        previousMonth=thisMonth-1;
        if previousMonth==0 % if this Month is January
            previousMonth=12; % previous month is December
            previousYear=thisYear-1; % year of previous month
        else
            previousYear=thisYear;  % year of previous month
        end
        
        
        
        
        % find month after this month
        nextMonth=thisMonth+1;
        if nextMonth==13 % if this Month is December
            nextMonth=1; % next month is January
            nextYear=thisYear+1; % year of previous month
        else
            nextYear=thisYear;  % year of next month
        end
%         
%         % find 2 months after this month
%         twomonthsAfter=thisMonth+2;
%         if twomonthsAfter==13 % if this Month is November
%             twomonthsAfter=1; % in 2 months, it's January
%             yearInTwoMonths=thisYear+1; % year in 2 months
%         elseif twomonthsAfter==14 % if this Month is December
%             twomonthsAfter=2; %  month is feb
%             yearInTwoMonths=thisYear+1; % year in two months
%         else
%             yearInTwoMonths=thisYear;  % year in two months
%         end
        
        
        % if there is a data entry for previous year and month:
        if any(YearMonthVolumePressureThisWell(:,1)==previousYear & YearMonthVolumePressureThisWell(:,2)==previousMonth)
            %  do nothing because that data entry exists
            %             volumeOfPreviousMonth=YearMonthVolumePressureThisWell(YearMonthVolumePressureThisWell(:,1)==previousYear & YearMonthVolumePressureThisWell(:,2)==previousMonth,3);
            
            % if injection rate didn't change from last month:
            if YearMonthVolumePressureThisWell(YearMonthVolumePressureThisWell(:,1)==previousYear & YearMonthVolumePressureThisWell(:,2)==previousMonth,3)==...
                    YearMonthVolumePressureThisWell(YearMonthVolumePressureThisWell(:,1)==thisYear & YearMonthVolumePressureThisWell(:,2)==thisMonth,3)
                % then last month's injection equals this months injection, so
                % move that time data point from the end of last month to the
                % end of this month.% this makes the datapoints not duplicate
                % when wells injct at a constant rate.
                datenumBarrelsPerDayThisWell(end,1)=[endDateThisMonthdatenum]; % move last months datenumber to this months datenumber, rate doesn't change
                continue % go to next month in for loop, just moved time point from previous month. 
            end
        else % there is no entry for previous month, so make last second of previous month data zero
            %             endDatePreviousMonthdecyear=decyear(previousYear,previousMonth,eomday(previousYear,previousMonth),23,59,59); % last second of previous month
            %             decyearBarrelsPerDayThisWell(end+1,[1:2])=[endDatePreviousMonthdecyear,0]; % make entry of 0 barrels previous month
            
            endDatePreviousMonthdatenum=datenum(previousYear,previousMonth,eomday(previousYear,previousMonth),23,59,59); % last second of previous month
            datenumBarrelsPerDayThisWell(end+1,[1:2])=[endDatePreviousMonthdatenum,0]; % make entry of 0 barrels previous month
        end
        
        %         decyearBarrelsPerDayThisWell(end+1,[1:2])=[startDateThisMonthdecyear,barrelsPerDayOverThisMonth]; % make entry at beginning of month  (green dot)
        %         decyearBarrelsPerDayThisWell(end+1,[1:2])=[endDateThisMonthdecyear,barrelsPerDayOverThisMonth]; % make entry at end of month  (green dot)
        
        datenumBarrelsPerDayThisWell(end+1,[1:2])=[endDateThisMonthdatenum,barrelsPerDayOverThisMonth]; % make entry at end of month  (black dot)
        
%         
%         if any(YearMonthVolumePressureThisWell(:,1)==nextYear & YearMonthVolumePressureThisWell(:,2)==nextMonth) % if there is a data entry for next month
%             % do nothing because that data entry exists
%         else % no entry for next month starting with 0
%             startDateNextMonthdecyear=decyear(nextYear,nextMonth,1,0,0,0); % beginning of first day of next month
%             %             decyearBarrelsPerDayThisWell(end+1,[1:2])=[startDateNextMonthdecyear,0]; % make entry of 0 barrels next month (green dot)
%             
%             if ~any(YearMonthVolumePressureThisWell(:,1)==yearInTwoMonths & YearMonthVolumePressureThisWell(:,2)==twomonthsAfter) % if there is no data entry for two months after
%                 startDateNextMonthdatenum=datenum(nextYear,nextMonth,1,0,0,0); % beginning of first day of next month
%                 datenumBarrelsPerDayThisWell(end+1,[1:2])=[startDateNextMonthdatenum,0]; % make entry at beginning of month  (black x)
%             end
%             
%         end
        
        
    end % done cycling over months of this well
    
    
    
    assumeInjectionExtendsToEndOfModel=true;
    
    if assumeInjectionExtendsToEndOfModel
        if nargin ==0 % for coding
            YearEndOfModel=2045;
        else
            YearEndOfModel=get(hDV.hdsldr(1),'max');% end of scaleruler.
        end
        %         maxyear=2045; % need a way to find this from max(slider years)
        %         decyearBarrelsPerDayThisWell(end,2)=[decyearBarrelsPerDayThisWell(end-1,2)]; % last second of last year in model
        %         decyearBarrelsPerDayThisWell(end+1,1:2)=[decyear(maxyear-1,12,31,23,59,59),decyearBarrelsPerDayThisWell(end,2)]; % last second of last year in model
        
%                 datenumBarrelsPerDayThisWell(end,1)=[endDateThisMonthdatenum]; % move last months datenumber to this months datenumber, rate doesn't change
% 

        datenumBarrelsPerDayThisWell(end,1)=[datenum(YearEndOfModel,1,0,0,0,0)];% last year in model

        
        
    end
    
    %     hDV.data.realWellData.decyearBarrelsPerDay{CycleWellsK,1}=decyearBarrelsPerDayThisWell; % save decyear and bblsPerDay that will be plotted
    hDV.data.realWellData.datenumBarrelsPerDay{CycleWellsK,1}=datenumBarrelsPerDayThisWell; % save datenumber and bblsPerDay , you need to do stairs(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2)) to plot this
%     datenumBarrelsPerDayThisWell
    
end % done cycling over all wells



% Clear any data of previous wells beyond hDV.data.nwells:
hDV.data.realWellData.APIs=hDV.data.realWellData.APIs(1:hDV.data.nwells,1);
hDV.data.realWellData.XEasting=hDV.data.realWellData.XEasting(1:hDV.data.nwells,1);
hDV.data.realWellData.YNorthing=hDV.data.realWellData.YNorthing(1:hDV.data.nwells,1);
hDV.data.realWellData.wellNames=hDV.data.realWellData.wellNames(1:hDV.data.nwells,1);
hDV.data.realWellData.wellTypes=hDV.data.realWellData.wellTypes(1:hDV.data.nwells,1);
hDV.data.realWellData.YearMonthVolumePressure=hDV.data.realWellData.YearMonthVolumePressure(1:hDV.data.nwells,1);
% hDV.data.realWellData.decyearBarrelsPerDay=hDV.data.realWellData.decyearBarrelsPerDay(1:hDV.data.nwells,1);
hDV.data.realWellData.datenumBarrelsPerDay=hDV.data.realWellData.datenumBarrelsPerDay(1:hDV.data.nwells,1);


% d=hDV.data.realWellData.decyearBarrelsPerDay




if nargin==0
    % check by plotting 
    % plot
    figure
    hold on
for dd=1:length(hDV.data.realWellData.datenumBarrelsPerDay)
  datenumBarrelsPerDayThisWell=hDV.data.realWellData.datenumBarrelsPerDay{dd,1};
[tt,QQ]=stairs(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2));
plot(datenumBarrelsPerDayThisWell(:,1),datenumBarrelsPerDayThisWell(:,2),'rx') % data points as red xs

% do this adjustment because value represents rate up to that point, not
% rate after that point as stair assumes. 
% see plottingInjectionStairData.m
QQ(1:2)=[];
QQ(end+1:end+2)=QQ(end);
plot(tt,QQ)
end

datetick('x')

    
end % end SpreadsheetStrings2WellData function











