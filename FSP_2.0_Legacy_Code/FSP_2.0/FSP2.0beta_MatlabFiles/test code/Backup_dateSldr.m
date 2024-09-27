% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%GUI Slider



function [sldr_handle ,sldr_txt]=dateSldr(hfig,PI,WM,HE,GW,hDV)
   
   % Extract input arguments
   cc.hfig     = hfig;
   cc.hCallObj = hDV;
   
   % Value range
   cc.minlim   = 2015;
   cc.maxlim   = 2045;
   dl = cc.maxlim - cc.minlim ; 
   cc.xleft    = cc.minlim;
   cc.xright   = cc.maxlim;
   cc.mindelta = 0.002;
   
   % Parameters for setting up sliders
   n     = 5; alpha = 1.2;
   
   % Time and slider button widths
   WD = (WM - (2 + 4 * n) * GW) / (2 + 3 * alpha);
   WS = alpha * WD;
   
   % Initial position and button sizes
   init = [(PI(1) + GW) PI(2) 0 0];
   szdb = [0 0 WD HE];
   szsb = [0 0 WS HE];
   
     % 'Year:' text
   post = init + szdb;
   uicontrol(cc.hfig,'style','text',                  ...
      'HandleVisibility','callback','Position',post,                   ...
      'HorizontalAlignment','left','string','Year:','fontsize',18, ...
      'fontunits','normalized','backgroundcolor',hDV.colors.white);
  
     % Display numerical year
   post = init + szdb + [WD/2 0 0 0];
   cc.htxtLsled = uicontrol(cc.hfig,'style','text',                  ...
      'HandleVisibility','callback','Position',post,'string',num2str(cc.minlim),   ...
      'HorizontalAlignment','center','fontsize',18,'fontunits','normalized',...
      'backgroundcolor',hDV.colors.white);
  

   % Slider
   post = init + [1.5*WD 0 0 0] + szsb;
   cc.hLsled = uicontrol(cc.hfig,'Style','slider','Tag','Lsled',       ...
      'String','Lsled','Min',cc.minlim,'Max',cc.maxlim,'Position',post,...
      'Value',cc.minlim,'HandleVisibility','callback', 'SliderStep',[1/dl 5/dl],  ...
      'BackgroundColor',[0.7 0.7 0.75],'Callback',{@sldrcall,cc.htxtLsled,hDV});
    
    %initialize
   sldr_handle=cc.hLsled;
   sldr_txt=cc.htxtLsled;
    
    function sldrcall(src,~,yrtxt,hDV)
        val = get(src,'value') ; val = floor(val) ;
        mival = get(src,'min') ; mxval = get(src,'max') ;
        if val<mival
            val=mival;
        elseif val>mxval
            val=mxval;
        end
        set(src,'value',val) ; 
        set(yrtxt,'string',num2str(val)) ; 
        
        %update all other sliders and text 
        set(hDV.hdsldr(:),'value',val) ; 
        set(hDV.hdsldr_txt(:),'string',num2str(val)) ;
        
        calcengine(hDV,hDV.currtab.name); %rerun calculations
        refreshplotdata(hDV,hDV.currtab.name); %plot data for selected panel refreshed
    end
  
   
   
end  % function dateSldr

