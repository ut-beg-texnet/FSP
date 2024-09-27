% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Standard data entry GUI
% modified by Rall for Probabilistic data entry
% and varying APhi/Mi or SH and Sh

function valsout = dataentryProbabilistic(hDV,title,txt,vals,noshow,deterministicValNums)

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
pos = [0.06 -0.07 .38 .05] ;
c2=[.53 0 -.2 0] ;% left box
c3=[.74 0 -.2 0] ; % right box
c4=[.34 0 -.2 0] ; % dropdown menus
valtxtH = zeros(length(vals),1) ;

strings1={'Plus/Minus'};
m=1 ;

boundlabel =uicontrol('parent',figureDataEntryProb,'style','text','string','Plus/Minus',...
    'pos',c2+pos+[0 .95-m*.08 0 0]+[.0 .035 0 -0.001],'fontsize',9,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'HorizontalAlignment','center','visible','on');

for k=1:1:length(txt)
    if ~any(k==noshow)
        
        uicontrol('parent',figureDataEntryProb,'style','text','string',txt{k},...
            'pos',pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1],'HorizontalAlignment','left');
        %  field
        valtxtH(k,1)= uicontrol('parent',figureDataEntryProb,'style','edit','string',num2str(vals(k,1)),...
            'pos',c2+pos+[0 .95-m*.08 0 0],'fontsize',12,'fontunits','normalized',...
            'backgroundcolor',[1 1 1],'visible','on' );
        m=m+1 ;
    else
        valtxtH(k,1)= uicontrol('parent',figureDataEntryProb,'style','edit','string',num2str(vals(k)),'visible','off',...
            'pos',[0 0 .1 .1]);
    end
end


% Add text box to specify stress model being used
% 
str_model_txt = 'Undefined stress model'; % This line should never be shown!
if hDV.data.stress.aphi.use==0
    str_model_txt = 'Gradient based stress model is being used                               ';
elseif hDV.data.stress.aphi.use==11
    str_model_txt = 'A-Phi stress model is being used                                        ';
elseif hDV.data.stress.aphi.use==12
    str_model_txt = 'Modified A-Phi stress model with min horiz stress gradient is being used';
end
uicontrol('parent',figureDataEntryProb,'style','text','string',str_model_txt,...
    'pos',[0.1 .92 0.6 0.04],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1],'HorizontalAlignment','left');


%OK button
uicontrol('parent',figureDataEntryProb,'style','pushbutton','string','OK',...
    'pos',[.25 .02 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcallProb,hDV,valtxtH,noshow,figureDataEntryProb,deterministicValNums});
valsout=0;



    function crf(~,~,figureDataEntryProb)
        hDV.tempvec=[];
        delete(figureDataEntryProb);
        modal(hfig,'on');
    end


%this must tie back to callbackdatabuts.m
    function butcallProb(~,~,hDV,htxt,noshow,figureDataEntryProb,deterministicValNums)
        
        %button color change to red
        resetprobpanel(hDV);
        
        vec=zeros(size(htxt));
        for n=1:1:length(htxt)
            p=1;% used to cycle over if multiple columns
            if strcmpi('on',get(valtxtH(n,p),'visible')) && ~any(n==noshow) % if field value shown
                vec(n,p)=str2double(get(htxt(n,p),'string')); % get field value
                
                % make sure no negative numbers
                if vec(n,p)<0  || isnan(vec(n,p))
                    edlgbox = errordlg(cat(2,'check ',get(htxt(n,p),'string'),' you entered a negative or non number plus/minus value: ',num2str(vec(n,p)),...
                        ', FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution'));
                    centerFigure(hDV.hfig,edlgbox);
                    return
                end
                if any(n==[1,2,3,5,10,12,13,14]) % these are the values that can't  be negative in any iteration
                    if vec(n,p) >= deterministicValNums(n)
                        edlgbox = errordlg(cat(2,'check ',get(htxt(n,p),'string'),' you entered a plus/minus value: ',num2str(vec(n,p)),...
                            ', that is greater than the deterministic value. FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution, but this parameter can''t be negative'));
                        centerFigure(hDV.hfig,edlgbox);
                        return
                    end
                end
                
                if n==13 && hDV.data.stress.aphi.use>0
                    if (vec(n,p) + deterministicValNums(n))>3
                        edlgbox = errordlg(cat(2,'check ',get(htxt(n,p),'string'),' you entered a plus/minus value: ',num2str(vec(n,p)),...
                            ', that when added to the deterministic aphi is greater than 3 This is impossible See Hurd and Zoback 2012 for definition of aPhi. '));
                        centerFigure(hDV.hfig,edlgbox);
                        
                        return
                    end
                end
            end
            
        end
        
        
        flag = checkProbabilisticEntries(hDV,vec,txt,noshow); % flag=1 if no issues
        
        if flag
            hDV.data.sigvals=vec(1:12); % store all plusminus values except A-Phi
            if  hDV.data.stress.aphi.use
                hDV.data.stress.aphi.sigvals(1)=vec(13); % store value for A-phi
            end
            delete(figureDataEntryProb);
            modal(hfig,'on');
        end
        
    end

% check probabilistic data distribution bound entries
    function [flag] = checkProbabilisticEntries(hDV,vec,txt,noshow)
        flag=true; % flag=1 if no issues
        % by Rall
        % still need to check relative stress magnitudes etc, and if stress state
        % is possible.
        for kj=1:length(vec);
            if vec(kj)<0 &&  ~any(kj ==noshow)
                edlgbox = errordlg(cat(2,'check ',txt{kj},' you entered a negative plus/minus value: ',num2str(vec(kj)),...
                    ', FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution'));
                centerFigure(hDV.hfig,edlgbox);
                flag=false;
                
            end
        end
    end

end




