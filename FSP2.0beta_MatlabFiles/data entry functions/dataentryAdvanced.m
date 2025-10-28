% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Standard data entry GUI

function valsout = dataentryAdvanced(hDV,title,txt,vals,i,noshow)

if nargin==5
    noshow=0; %entries which we dont want to show
end

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
f = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .4*hpos(3) .7*hpos(4)];
set(f,'position',pos) ;
set(f,'MenuBar','none');
set(f,'Name',title);
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'PaperPositionMode','auto');
set(f,'closerequestfcn',@crf) ;
centerFigure(hfig,f);
set(f,'visible','on') ; 
pos = [0.05 0 .45 .05] ;  c2=[.6 0 -.15 0] ;
valtxtH = zeros(length(vals),1) ;

toolTipStrings5={'Map domain limit';'Map domain limit';'Map domain limit';'Map domain limit';...
    'Dendity of injected fluid';'Gravitational constant g';'Fluid Viscosity';'Fluid Compressibility';'Injection Formation Compressibility';...
    'Zero doesn''t set a seed, positive integer makes the same "random" occurrence'};

m=1 ;
for k=1:1:length(txt)
    if ~any(k==noshow)
        uicontrol('parent',f,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.07 0 0],'fontsize',12,'fontunits','normalized',...
            'tooltipstring',toolTipStrings5{k},...
            'backgroundcolor',[1 1 1]);
        
        valtxtH(k)= uicontrol('parent',f,'style','edit','string',num2str(vals(k)),...
            'pos',c2+pos+[0 .95-m*.07 0 0],'fontsize',12,'fontunits','normalized',...
             'tooltipstring',toolTipStrings5{k},...
            'backgroundcolor',[1 1 1]);
        m=m+1 ;
    else
        valtxtH(k)= uicontrol('parent',f,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
end
%OK button
uicontrol('parent',f,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,i});

valsout=0;

%debug graphicsButton - added in V0.98, only reason this isn't default
%dataentry function anymore
if 1
    uicontrol('parent',f,'style','text','string','Choose Renderer:',...
        'pos',pos+[0 .95-m*.07 0 0],'fontsize',12,'fontunits','normalized',...
        'backgroundcolor',[1 1 1],'tooltipstring',...
        'some graphics renderers are incompatible with some grapics cards, if necessary, choose an alternative');
    
    uicontrol('parent',f,'Style','popup','String','openGL Software|openGL hardware|zBuffer|manual',...
        'pos',c2+pos+[0 .95-m*.07 0 0],'fontsize',12,'fontunits','normalized',...
        'callback',{@toggleRenderer,hDV},'enable','on','value',1, 'backgroundcolor',[1 1 1],...
        'tooltipstring','some graphics renderers are incompatible with some grapics cards, if necessary, choose an alternative');
end
             
    function crf(~,~)
        hDV.tempvec=[];
        delete(f);
        modal(hfig,'on');
    end

    %this must tie back to callbackdatabuts.m
    function butcall(~,~,hDV,htxt,i)
        vec=zeros(length(htxt),1);
        for n=1:1:length(htxt)
            vec(n)=str2double(get(htxt(n),'string'));
        end
        vec
        switch i
            case 1
                hDV.data.stress.vals=vec;
            case 2
                hDV.data.reservoir.vals=vec;
            case 3
                hDV.data.inject.vals=vec;
            case 4
                hDV.data.fault.vals=vec;
            case 5
                if ~all(hDV.data.adv.vals==vec)% then some data was changed
                    resetButtonsRed(hDV,[1,2,3,4])            % Set all buttons red
                end
                hDV.data.adv.vals=vec;
            case 6
                hDV.data.sigvals=vec;
                updateinputbar(hDV);
        end
        delete(f);
        modal(hfig,'on');
    end


end




