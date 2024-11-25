function callinthelp(~,~,hDV)

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
f = figure('color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .5*hpos(3) .5*hpos(4)];
set(f,'position',pos) ;
set(f,'MenuBar','none');
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'closerequestfcn',@crf) ;
centerFigure(hfig,f); 
imshow('helpint.png','border','tight');
set(f,'visible','on') ; 
    function crf(~,~)
        delete(gcf);
        modal(hfig,'on');
    end


end
