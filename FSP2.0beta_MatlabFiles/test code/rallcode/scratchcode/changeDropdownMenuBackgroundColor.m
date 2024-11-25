% make dorodownmnu color background scratsh 


figure
% make dropdownMenuOptions background colored
crramp='jet';

menuOptions={'aasd';'fdsab';'casdf';'dasdf';'fffffas'};
endc = eval([char(crramp) '(',num2str(length(menuOptions)),')']); % crramp is a string 
%defining the colormap
for a = 1:size(endc,1)
    pue{a} = ['<HTML><BODY bgcolor="rgb(' num2str(endc(a,1)*255) ',' ...
        num2str(endc(a,2)*255) ',' num2str(endc(a,3)*255) ')">  '...
        menuOptions{a},'    </BODY></HTML>'];
end
uicontrol('Units','normalized',...
    'Position',[.5,.1,.4,.5],'HorizontalAlignment','left',...
  'Style','popup','String',pue); % set to color boxes


