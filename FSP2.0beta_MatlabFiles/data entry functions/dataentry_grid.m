% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface


% Grid data entry GUI


function valsout = dataentry_grid(hDV,title,txt,vals,i)

hfig = hDV.hfig;
modal(hfig,'off');

[M,N] = size(vals);  %dimensions of data entry matrix; rows are individual wells
valtxtH = zeros(M,N) ;

hpos = get(hfig,'position') ;
f = figure('tag',title,'color',[1 1 1],'units','pixels');
pos = [0 0 .7*hpos(3) .9*hpos(4)];
set(f,'position',pos) ;
set(f,'MenuBar','none');
set(f,'Name',title);
set(f,'NumberTitle','off');
set(f,'DefaultUicontrolUnits','normalized');
set(f,'DefaultUicontrolFontsize',14);
set(f,'PaperPositionMode','auto');
set(f,'closerequestfcn',@crf) ;
centerFigure(hfig,f);

uicontrol('parent',f,'style','text','string','Number of wells:',...
    'pos',[.05 .9 .4 .05],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);

pos0=[.05 0 .1 .05];
for k=1:1:M
    uicontrol('parent',f,'style','text','string',char(65+k-1),...
    'pos',pos0+[0 .8-k*.08 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1]);
end

wd=.8/length(txt) ; pos = [0 .8 wd .05] ;  
for k=1:1:N
    uicontrol('parent',f,'style','text','string',txt{k},...
        'pos',pos+[.15+(k-1)*wd 0 0 0],'fontsize',12,'fontunits','normalized',...
        'backgroundcolor',[1 1 1]);
end



for m=1:1:M
    for k=1:1:N
        posoff = [0 -m*.08  0 0] ;
        valtxtH(m,k)= uicontrol('parent',f,'style','edit','string',num2str(vals(m,k)),...
            'pos',pos+[.15+(k-1)*wd 0 0 0]+posoff,'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1]);
    end
end


hNwells = uicontrol('parent',f,'style','edit','string',num2str(hDV.data.nwells),...
    'pos',[.45 .9 .2 .05],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'callback',{@nwellscall,valtxtH,M});

%setup to start
nwellscall(hNwells,[],valtxtH,M);


%OK button
uicontrol('parent',f,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,i});

valsout=0;

    function crf(~,~)
        hDV.tempvec=[];
        delete(gcf);
        modal(hfig,'on');
    end

    %this must tie back to callbackdatabuts.m
    function butcall(~,~,hDV,htxt,i)
        [P,Q]=size(htxt); vec=zeros(P,Q);
        
        for p=1:P
            for q=1:Q
                vec(p,q)=str2double(get(htxt(p,q),'string'));
            end
        end
        
        switch i
            case 3
                hDV.data.inject.vals=vec;
        end
        delete(gcf);
        modal(hfig,'on');
    end


    function nwellscall(src,~,hvalstxt,maxval)
        
        val = str2double(get(src,'string')) ;
        val=floor(val);
        if val<1
            val=1;
        elseif val>maxval
            val=maxval;
        end
        set(src,'string',num2str(val)); 
        
        hDV.data.nwells = str2double(get(src,'string')) ; 
        
        for j=1:hDV.data.nwells
            set(hvalstxt(j,:),'visible','on') ; 
        end
        
        for j=hDV.data.nwells+1:M
            set(hvalstxt(j,:),'visible','off') ; 
        end
    end


end




