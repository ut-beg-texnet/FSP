% color text of dropdown menu or listbox menu 
% By Rall Walsh
% Stanford
% inspired by
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/163357
% November 2016

% need to improve with number of faults check, if FSP has a value check,
% etc


% cv=hDV.plotdata.pint.fsp; % get fault FSPs



function  probHydDropdownTxtClr(handleOfMenu,hDV,parameter2Colorby,minint,maxint,colormapToUse,addFspBinary)
if nargin==0
    hDV=evalin('base', 'hSV'); % get hDV to test?
    parameter2Colorby=hDV.plotdata.pint.fsp;
    handleOfMenu=hDV.plotdata.ListboxFaultSelector;
    minint=hDV.plotdata.minint;
    maxint=hDV.plotdata.maxint;
    colormapToUse=hDV.cmapGYR;
    addFspBinary=1;
elseif nargin==2
    parameter2Colorby=hDV.plotdata.pint.fsp;
    minint=hDV.plotdata.minint;
    maxint=hDV.plotdata.maxint;
    colormapToUse=hDV.cmapGYR;
    addFspBinary=1;
end


nfaults=hDV.data.fault.vals(1); % get fault count
if length(parameter2Colorby)~=nfaults
    msgWindow1=msgbox('error coloring text in probHydDropdownTxtClr, fault FSP vector length is different from nfaults.', ' Internal code error','error'); % this should never be executed
    centerFigure(hDV.hfig,msgWindow1);
end

% faultFSP=hDV.plotdata.pint.fsp;
% colormap=hDV.cmapGYR; % get trafficlight colormap
endc=zeros(nfaults,3); % preallocate
% endc=colormap(1:nfaults,:);  % rgb of each fault
% cv=hDV.plotdata.pint.fsp; % get fault FSPs

if isnan(minint) || isnan(maxint)
    minint=0 ; maxint=1;
end

xstr = cell(nfaults+1,1) ; xstr{1}='All Faults' ;
for k=1:1:nfaults
    fspToColor=parameter2Colorby(k); % value
    cl = getcolor(colormapToUse,fspToColor,minint,maxint) ; % find RGB color
    endc(k,1:3)=cl;
    if addFspBinary
        xstr{k+1} = ['<HTML><BODY bgcolor="rgb(' num2str(endc(k,1)*255) ',' ...
        num2str(endc(k,2)*255) ',' num2str(endc(k,3)*255) ')">  '...
        'Fault #',num2str(k),', ',num2str(fspToColor,'%0.2f'),' FSP    </BODY></HTML>'];
    else
         xstr{k+1} = ['<HTML><BODY bgcolor="rgb(' num2str(endc(k,1)*255) ',' ...
        num2str(endc(k,2)*255) ',' num2str(endc(k,3)*255) ')">  '...
        'Fault #',num2str(k),'    </BODY></HTML>'];
  
    end
end
set(handleOfMenu,'string',xstr) ;


end



















