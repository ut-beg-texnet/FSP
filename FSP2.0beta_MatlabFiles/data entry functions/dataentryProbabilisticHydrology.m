% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Standard data entry GUI
% modified by Rall for Probabilistic data entry
% and varying APhi/Mi or SH and Sh

function valsout = dataentryProbabilisticHydrology(hDV,title,txt,vals,noshow)

if nargin==4
    noshow=0; %entries which we dont want to show
end

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
figureDataEntryProb = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .6*hpos(3) .7*hpos(4)];
set(figureDataEntryProb,'position',pos) ;
set(figureDataEntryProb,'MenuBar','none');
set(figureDataEntryProb,'Name',title);
set(figureDataEntryProb,'NumberTitle','off');
set(figureDataEntryProb,'DefaultUicontrolUnits','normalized');
set(figureDataEntryProb,'DefaultUicontrolFontsize',14);
set(figureDataEntryProb,'PaperPositionMode','auto');
set(figureDataEntryProb,'closerequestfcn',{@crf,figureDataEntryProb}) ;
centerFigure(hfig,figureDataEntryProb);
set(figureDataEntryProb,'visible','on') ;
pos = [0.06 0 .48 .05] ; % pos = [0.06 0 .38 .05] ;
c2=[.53 0 -.2 0] ;% left box
c3=[.74 0 -.2 0] ; % right box
c4=[.34 0 -.2 0] ; % dropdown menus
valtxtH = zeros(length(vals),1) ;


% radio buttons
radioButtonGroup=uibuttongroup('Position', [.3 .9 .4 .1],'backgroundcolor',[1 1 1],...
      'parent',figureDataEntryProb,'fontsize',hDV.ftsz,...
      'BorderType','none','fontunits','normalized');

radioButton1=uicontrol(radioButtonGroup,'Style','radiobutton',...
    'Position', [0.01 .45 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Probabilistic Hydrology','fontunits','normalized',...
    'ToolTipString','Select run probabilistic hydrology, which means making lots of hydrologic models');
  
radioButton2=uicontrol(radioButtonGroup,'Style','radiobutton',...
    'Position', [0.01 .0 .98 .5],'backgroundcolor',[1 1 1],...
    'String','Deterministic Hydrology','fontunits','normalized',...
    'ToolTipString','Select to just make one preferred hydrologic model');

probabilistic1VsDeterministic2= decideDeterministicOrProbHydrology(hDV);

    switch probabilistic1VsDeterministic2
        case 2 % deterministic 
             set(radioButtonGroup,'SelectedObject',radioButton2) 
             evt.NewValue=radioButton2;
        otherwise  % probabilistic
            set(radioButtonGroup,'SelectedObject',radioButton1)
            evt.NewValue=radioButton1;
    end

% type faults table
enterSigmasPanel = uipanel('parent',figureDataEntryProb,'backgroundcolor','white',    ...
    'visible','on','Position',[0 0.15 1 .75],'BorderType','none');

strings1={'Plus/Minus'};
m=1 ;

boundlabel =uicontrol('parent',enterSigmasPanel,'style','text','string','Plus/Minus:',...
    'pos',c2+pos+[0 .95-m*.08 0 0]+[.0 .035 0 -0.001],'fontsize',9,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'HorizontalAlignment','center','visible','on');

boundlabel2 =uicontrol('parent',enterSigmasPanel,'style','text','string','Change Computations?',...
    'pos',[ 0.5700,0.300, 0.41, 0.0490],'fontsize',9,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'HorizontalAlignment','center','visible','on');


for k=1:1:length(txt)
    if ~any(k==noshow)
        
        uicontrol('parent',enterSigmasPanel,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1],'HorizontalAlignment','left');
        %  field
        valtxtH(k,1)= uicontrol('parent',enterSigmasPanel,'style','edit','string',num2str(vals(k,1)),...
            'pos',c2+pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1],'visible','on','enable','on');
        m=m+1 ;
    else
        valtxtH(k,1)= uicontrol('parent',enterSigmasPanel,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
    if k==7 ;m=m+1; end % space bootstrap edit
end

set(radioButtonGroup,'SelectionChangeFcn', {@radioMenuProbHydToggle,hDV,radioButton1,radioButton2,valtxtH})

radioMenuProbHydToggle([],evt,hDV,radioButton1,radioButton2,valtxtH)

%OK button
uicontrol('parent',figureDataEntryProb,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcallProb,hDV,valtxtH,noshow,figureDataEntryProb,txt});
valsout=0;

    function crf(~,~,figureDataEntryProb)
        hDV.tempvec=[];
        delete(figureDataEntryProb);
        modal(hfig,'on');
    end

%this must tie back to callbackdatabuts.m
    function butcallProb(~,~,hDV,htxt,noshow,figureDataEntryProb,txt)

        %button color change to red
        resetprobpanel(hDV);
        
        vec=zeros(size(htxt));
        for n=1:1:length(htxt)
            p=1;% used to cycle over if multiple columns
            if strcmpi('on',get(valtxtH(n,p),'visible')) && ~any(n==noshow) % if field value shown
                vec(n,p)=str2double(get(htxt(n,p),'string')); % get field value
            end
            
        end
        
        flag =   checkProbabilisticEntriesHyd(hDV,noshow,vec,txt); % flag=1 if no issues
        if flag
            set(hDV.plotdata.hydprob.hbutMC,'backgroundcolor','r')
            set(hDV.plotdata.pprob.summaryPlotButton,'backgroundcolor',hDV.colors.red)
            hDV.data.probHydrology.sigvals=vec; % store
            delete(figureDataEntryProb);
            modal(hfig,'on');
        end
    end

end



function radioMenuProbHydToggle(~,evt,hDV,radioButton1,radioButton2,valtxtH)

switch evt.NewValue
    
    case radioButton1 % probabilistic
        set(valtxtH(:),'enable','on')
        hDV.data.probHydrology.probabilistic1VsDeterministic2=1;
    case radioButton2 % deterministic
        set(valtxtH(:),'enable','off')
        hDV.data.probHydrology.probabilistic1VsDeterministic2=2;
        
        if hDV.data.reservoir.importHydrology   % imported modflow model? 1=yes
            set(radioButton1,'enable','off')
            set(radioButton2,'ToolTipString','Because You imported a deterministic hydrologic model, you can''t run probabilistic hydrology')
        end
        
        
end



end












