% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%GUI Slider
% Modified by Rall Walsh to make year box editable too, Stanford 2016


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





% Slider
post = init + [1.5*WD 0 0 0] + szsb;
cc.hLsled = uicontrol(cc.hfig,'Style','slider','Tag','Lsled',       ...
    'String','Lsled','Min',cc.minlim,'Max',cc.maxlim,'Position',post,...
    'Value',cc.minlim,'HandleVisibility','callback', 'SliderStep',[1/dl 5/dl],  ...
    'BackgroundColor',[0.7 0.7 0.75],'Callback',{@sldrcall,hDV});
%     'BackgroundColor',[0.7 0.7 0.75],'Callback',{@sldrcall,cc.htxtLsled,hDV});

% Display numerical year
post = init + szdb + [WD/2 0 0 0];
cc.htxtLsled = uicontrol(cc.hfig,'style','edit',                  ...
    'HandleVisibility','callback','Position',post,'string',num2str(cc.minlim),   ...
    'HorizontalAlignment','center','fontsize',18,'fontunits','normalized',...
    'backgroundcolor',hDV.colors.white,...
    'Callback',{@yearEditCall,cc,hDV});

%initialize
sldr_handle=cc.hLsled;
sldr_txt=cc.htxtLsled;

      function yearEditCall(src,~,cc,hDV) % callback if year string is edited - call slider callback
        
        valstr = get(src,'string')  ;
        val = floor(str2double(valstr)) ;
        mival = get(hDV.hdsldr(1),'min');
        mxval = get(hDV.hdsldr(1),'max');
        if val<mival
            val=mival;
        elseif val>mxval
            val=mxval;
        end
        set(hDV.hdsldr(:),'value',val)
        sldrcall(cc.hLsled,src,hDV) % slider callback
        
    end


end  % function dateSldr

