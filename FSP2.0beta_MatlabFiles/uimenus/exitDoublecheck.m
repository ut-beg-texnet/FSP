% do you really want to quit or overwrite?


function exitDoublecheck(src,~,hDV)
hfig = hDV.hfig;
modal(hfig,'off');

title=['Are you sure you want to ',get(src,'tag'),'?'];
hpos = get(hfig,'position') ;
f = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .4*hpos(3) .1*hpos(4)];
set(f,'position',pos) ;
set(f,'MenuBar','none');
set(f,'Name',title);
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'PaperPositionMode','auto');
set(f,'closerequestfcn',{@crf,f,hDV}) ;
centerFigure(hfig,f);
set(f,'visible','on') ;

uicontrol('parent',f,'style','text','string',title,...
            'pos',[.05 .5 .9 .45],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
        
uicontrol('parent',f,'style','pushbutton','string','Yes!',...
    'pos',[.1 .1 .4 .4],'fontsize',12,'fontunits','normalized',...
    'tag',get(src,'tag'),'callback',{@yesIDo,f,hDV});

uicontrol('parent',f,'style','pushbutton','string','No, Don''t',...
    'pos',[.5 .1 .4 .4],'fontsize',12,'fontunits','normalized',...
    'tag',get(src,'tag'),'callback',{@crf,f,hDV});


end


% close request function: if cancel
   function crf(~,~,f,hDV)
        delete(f);
        modal(hDV.hfig,'on');
   end
   
   % if yes is pushed
   function yesIDo(src,~,f,hDV)
   switch get(src,'tag') % get tag of original button, which is text
       case 'quit'
           delete(hDV) % exit program
       case 'overwrite your project'
           modal(hDV.hfig,'on')
           callbackloaddata([],[],hDV) % load default data
   end
   delete(f)
   
   end