function [strdate,ret] = ymdhms2date(str,~)

%--------------------------------------------------------------------------
% (c) ExxonMobil Upstream Research Company, 2003+
%     Drilling & Subsurface Technology Division
%     Well Construction Section
%  Authors: Bailey, JR
%
%  Extract from input string the Matlab date/time code
%--------------------------------------------------------------------------

strdate = [];
ret = -1;

sl = strfind(str,'\');  % sl = findstr(str,'\');  % Payette
if ( ~isempty(sl) )
   [~,ff] = fileparts(str);
   dash = dirute('getdash',ff);   %findstr(str,'-');
else
   dash = strfind(str,'-');  % dash = findstr(str,'-');  % Payette
   ff = str;
end

if ( length(dash) < 2 || isempty(dash) )
   return;
end

duh = dash(end-1);

if ( dash(end) - duh == 3 )
   
   yyyy = str2double(ff(duh-4:duh-1));
   mo = str2double(ff(duh+1:duh+2));
   dd = str2double(ff(duh+4:duh+5));
   
   hh = [];  mm = [];  ss = [];
   sp = dash(end)+3;
   if ( length(ff) >= sp )
      if ( char(ff(sp:sp)) == ' ' )
         bal = ff(sp+1:end);
         try
            if ( length(bal) < 8 )
               if ( length(bal) == 2 )
                  hh = str2double(bal(1:2));
               elseif ( length(bal) == 5 )
                  hh = str2double(bal(1:2));
                  mm = str2double(bal(4:5));
               end
            elseif ( bal(3:3) == '.' || bal(3:3) == ':' )
               hh = str2double(bal(1:2));
               if ( bal(6:6) == '.' || bal(6:6) == ':' )
                  mm = str2double(bal(4:5));
                  ss = bal(7:8);
                  if ( ~strcmpi(ss,'st') ...
                        && ~strcmpi(ss,'ds') )   % block out sta,dsf files
                     ss = str2double(ss);
                  else
                     ss = 0;
                  end
               end
            end
         catch ME
            msg = ME.message;
            disp(['ymdhms2date(55)[',str,']::',msg]);
         end
      end
   end
   
   if ( isempty(hh) ),  hh = 0;  end
   if ( isempty(mm) ),  mm = 0;  end
   if ( isempty(ss) ),  ss = 0;  end
   strdate = datenum(yyyy,mo,dd,hh,mm,ss);
   
   % need to round strdate to nearest ms
   
   msperday = 1000*60*60*24;
   strdate = round(strdate*msperday)/msperday;
   
   ret = 0;
   
end
