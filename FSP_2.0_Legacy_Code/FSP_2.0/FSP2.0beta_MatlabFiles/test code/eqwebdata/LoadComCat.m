function [year, month, day, hour, minute, sec, lat, long, depth, mag, magType] = LoadComCat(startTime,endTime,minMagnitude,varargin)
%LOADCOMCAT         Batch query ComCat searches to get around NEIC 20,000 event limit
%        [YEAR, MONTH, DAY, HOUR, MINUTE, SEC, LAT, LONG, DEPTH,
%        MAG, MAGTYPE] = LOADCOMCAT(STARTTIME,ENDTIME,MINMAGNITUDE)
%        returns results of a global catalog search within the time
%        frame STARTTIME and ENDTIME, for events with magnitude
%        greater than MINMAGNITUDE.
%
%        [YEAR, MONTH, DAY, HOUR, MINUTE, SEC, LAT, LONG, DEPTH,
%        MAG, MAGTYPE] = LOADCOMCAT(STARTTIME,ENDTIME,MINMAGNITUDE,
%        [MINLAT MAXLAT MINLON MAXLON]) performs a search within the
%        specified lat/long box.
%
%        STARTTIME and ENDTIME must be entered in serial date
%        number format, e.g. STARTTIME = datenum('2014-01-01 00:00:00')
%
%        For searches that will return more than 20,000 events,
%        this code will perform multiple ComCat searches, in
%        series, and return the aggregated results.  Results will
%        be sorted in time, from oldest to newest events.
%
%        Uses the ComCat search API.  For more info see: 
%        http://comcat.cr.usgs.gov/fdsnws/event/1/
%
%        Authors: Morgan Page and Justin Rubinstein
%                 U. S. Geological Survey
%        Last modified: March 2015
  

if(isempty(varargin))
    minlat=-90;maxlat=90;minlon=-180;maxlon=180;
else
    jj=varargin{1};
    minlat=jj(1);maxlat=jj(2);minlon=jj(3);maxlon=jj(4);
end
  

% URL for ComCat "count" method
url='http://earthquake.usgs.gov/fdsnws/event/1/count';
clear params;
params{1}='starttime'; params{3}= 'endtime'; 
params{5}='minmagnitude';  params{6}= num2str(minMagnitude);
params{7}='minlongitude';params{8}=num2str(minlon);
params{9}='maxlongitude';params{10}=num2str(maxlon);
params{11}='minlatitude';params{12}=num2str(minlat);
params{13}='maxlatitude';params{14}=num2str(maxlat);
params{15}='eventtype';params{16}='earthquake'; % exclude mine blasts!
startTimes=startTime; endTimes=endTime;


% Find number of events in search criteria
clear count;
% keyboard
for i=1:length(startTimes)
  params{2}=datestr(startTimes(i)); 
  params{4}=datestr(endTimes(i));
  count(i) = str2num(urlread(url,'get',params));
  display(['Total catalog = ' num2str(count) ' events'])
end


% While some searches go over 20,000 event limit, continue to
% divide time periods in half
while max(count)>20000
  indicesOverLimit=find(count>20000);
  for ind=indicesOverLimit
    endTimes(end+1)=endTimes(ind);
    startTimes(end+1)=startTimes(ind)+(endTimes(ind)-startTimes(ind))/2;
    endTimes(ind)=startTimes(ind)+(endTimes(ind)-startTimes(ind))/2;
  end
  
  clear count;
  for i=1:length(startTimes)
    params{2}=datestr(startTimes(i)); 
    params{4}=datestr(endTimes(i));
    count(i) = str2num(urlread(url,'get',params));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have a set of time periods that will keep each query below 20,000 events
% Now batch load catalogs!

% URL for ComCat "query" method
url='http://earthquake.usgs.gov/fdsnws/event/1/query';
params{end+1}='format'; params{end+1}='csv';
params{end+1}='orderby'; params{end+1}='time-asc';
  
startTimes=sort(startTimes); endTimes=sort(endTimes);
data=[];
for i=1:length(startTimes)
  params{2}=datestr(startTimes(i));
  params{4}=datestr(endTimes(i));
  newdata=urlread(url,'get',params);
  newdata=newdata(min(find(newdata==char(10)))+1:end); % Throw out first line (header)
  data = [data newdata];  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Catalog is loaded in "data" vector
% Format:
% time,latitude,longitude,depth,mag,magType,nst,gap,dmin,rms,net,id,updated,"place",type
% Now let's parse the data into some MatLab variables
% Note the "place" field can have commas and spaces within it, so
% we can't easily use textscan()


% Get rid of spaces so we can use strread()
d=data;
d(find(d==' '))='_';
% Break up into individual event lines
[eventLine]=strread(d,'%s\n');


% Initialize variables
year=zeros(1,length(eventLine));
month=zeros(1,length(eventLine));
day=zeros(1,length(eventLine));
hour=zeros(1,length(eventLine));
minute=zeros(1,length(eventLine));
sec=zeros(1,length(eventLine));
lat=zeros(1,length(eventLine));
long=zeros(1,length(eventLine));
depth=zeros(1,length(eventLine));
mag=zeros(1,length(eventLine));
magType=cell(1,length(eventLine));


% Loop through events, save info!
for i=1:length(eventLine)
  eventString=eventLine{i};
  year(i)=str2num(eventString(1:4));
  month(i)=str2num(eventString(6:7));
  day(i)=str2num(eventString(9:10));
  hour(i)=str2num(eventString(12:13));
  minute(i)=str2num(eventString(15:16));
  sec(i)=str2num(eventString(18:22));  
  commaLocations=find(eventString==',');
  lat(i)=str2num(eventString(commaLocations(1)+1:commaLocations(2)-1));
  long(i)=str2num(eventString(commaLocations(2)+1:commaLocations(3)-1));
  
  if(isempty(str2num(eventString(commaLocations(3)+1:commaLocations(4)-1))))
      depth(i)=NaN;
  else
      depth(i)=str2num(eventString(commaLocations(3)+1:commaLocations(4)-1));
  end

  if isempty(str2num(eventString(commaLocations(4)+1:commaLocations(5)-1)))
    mag(i)=NaN; % Special case for some old events that have a blank magnitude field
  else
    mag(i)=str2num(eventString(commaLocations(4)+1:commaLocations(5)-1)); 
  end
  magType{i}=eventString(commaLocations(5)+1:commaLocations(6)-1); 
end







  
  
  
  
  
  
  
  
  
  
