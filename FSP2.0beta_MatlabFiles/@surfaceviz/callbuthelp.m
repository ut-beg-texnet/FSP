function callbuthelp(src,~,hDV,fileTitle,widthHeightProportions)

if nargin ==3
    fileTitle='help.jpg';
    widthHeightProportions=[0.5 , 0.5 ]; % fraction of width and height of main figure that help figure will be
elseif nargin == 4
    widthHeightProportions=[0.5 , 0.5 ];
end
hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;

f = figure('color',[1 1 1],'units','pixels','visible','off');
pos = [hpos(1), hpos(2), round(widthHeightProportions(1)*hpos(3)), round( widthHeightProportions(2)*hpos(4))];

set(f,'MenuBar','none');
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'closerequestfcn',{@crf,hfig,f}) ;
imshow(fileTitle,'border','tight');
set(f,'position',pos) ;
set(f,'visible','on') ; 


    function crf(~,~,hfig,f)
        delete(f);
        modal(hfig,'on');
    end


end
