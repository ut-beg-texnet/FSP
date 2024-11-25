% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Standard data entry GUI
% modified by Rall for Probabilistic data entry
% and varying APhi/Mi or SH and Sh

function valsout = dataentryProbabilistic(hDV,title,txt,vals,noshow)

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
pos = [0.06 0 .38 .05] ;
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

%OK button
uicontrol('parent',figureDataEntryProb,'style','pushbutton','string','OK',...
    'pos',[.25 .05 .5 .1],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcallProb,hDV,valtxtH,noshow,figureDataEntryProb});
valsout=0;

    function crf(~,~,figureDataEntryProb)
        hDV.tempvec=[];
        delete(figureDataEntryProb);
        modal(hfig,'on');
    end

%this must tie back to callbackdatabuts.m
    function butcallProb(~,~,hDV,htxt,noshow,figureDataEntryProb)

        %button color change to red
        resetprobpanel(hDV);
        
        vec=zeros(size(htxt));
        for n=1:1:length(htxt)
            p=1;% used to cycle over if multiple columns
            if strcmpi('on',get(valtxtH(n,p),'visible')) && ~any(n==noshow) % if field value shown
                vec(n,p)=str2double(get(htxt(n,p),'string')); % get field value
            end
            
        end
        
        % edit vec if entered aphi and mu
        if  hDV.data.stress.aphi.use % if field value shown is Aphi or mu entry, calculate corresponding SHmax and Shmin

            dpth=hDV.data.stress.vals(5); % 
            % low value, deterministic value, high value
            SvDetGrad=hDV.data.stress.vals(1); %
            SvPlusMinusValue=vec(1);
            SvGrads=[SvDetGrad-SvPlusMinusValue,SvDetGrad,SvDetGrad+SvPlusMinusValue] ;
            SvVals=SvGrads.*dpth;
            PpDetGrad= hDV.data.stress.vals(6); %
            PpPlusMinusValue= vec(5); %
            PpGrads=[PpDetGrad-PpPlusMinusValue,PpDetGrad,PpDetGrad+PpPlusMinusValue] ;
            PpVals=PpGrads.*dpth;
            APhiDetVal =hDV.data.stress.aphi.vals(1); %
            APhiPlusMinusVal=vec(2);
            aPhiVals=[APhiDetVal-APhiPlusMinusVal,APhiDetVal,APhiDetVal+APhiPlusMinusVal] ;
            muDetVal = hDV.data.stress.aphi.vals(2); % 
            muPlusMinusVal= vec(3);
            muVals=[muDetVal-muPlusMinusVal,muDetVal,muDetVal+muPlusMinusVal] ;
            if any(muVals)<0.0001
                muVals(muVals<0.0001)=0.0001;
                msgWindow1=msgbox({cat(2,'Warning, your Mu uncertainty: ',num2str(muPlusMinusVal),' is large compared to your mu value: ',num2str(muDetVal));...
                    cat(2,'so because mu can''t be 0 or negative, making lowest mu bound ',num2str(0.0001),'  ')},...
                    'Friction Bound impossibility warning','warn');
                centerFigure(hDV.hfig,msgWindow1);
            end
            
            if any(aPhiVals<0)
                aPhiVals(aPhiVals<0)=0;
                msgWindow2=msgbox({cat(2,'Warning, your APhi uncertainty: ',num2str(APhiPlusMinusVal),' is larger than your APhi value: ',num2str(APhiDetVal));...
                    cat(2,'so because Aphi can''t be negative, making lowest APhi bound ',num2str(0),'  ')},...
                    'APhi Bound impossibility warning','warn');
                centerFigure(hDV.hfig,msgWindow2);
            end
            if any(aPhiVals>3)
                aPhiVals(aPhiVals>3)=3;
                msgWindow3=msgbox({cat(2,'Warning, your APhi uncertainty: ',num2str(APhiPlusMinusVal),' plus your APhi value: ',num2str(APhiDetVal));...
                    cat(2,' allows values higher than 3. This is impossible(Hurd and Zoback 2012 Tectonophysics), so making highest APhi bound ',num2str(3),'  ')},...
                    'APhi Bound impossibility warning','warn');
                centerFigure(hDV.hfig,msgWindow3);
            end
                        
            
            matrixOfVals=[aPhiVals;SvVals;muVals;PpVals]; % rows are APhi,Sv,mu,Pp, Columns are low, mid, high
            AllCombinationsOfVals=combvec(aPhiVals,SvVals,muVals,PpVals);
            
            
            % preallocate
            AllSHMaxVals=zeros(size(AllCombinationsOfVals,2),1);
            AllShmaxVals=AllSHMaxVals; 
            
            % try each 3 values of a variable with each 3 of the other
            % variables, then take lowest and highest SH and Sh values as
            % bounds. 
            for CycleCases559=1:size(AllCombinationsOfVals,2)
              
                [AllSHMaxVals(CycleCases559,1),AllShmaxVals(CycleCases559,1)] = getHorFromAPhi(AllCombinationsOfVals(1,CycleCases559),AllCombinationsOfVals(2,CycleCases559),...
                    AllCombinationsOfVals(3,CycleCases559),AllCombinationsOfVals(4,CycleCases559),hDV);
                
            end
            
            % calculate deterministic value
            [middleDeterministicSHmax,middleDeterministicShmin] = getHorFromAPhi(aPhiVals(2),SvVals(2),muVals(2),PpVals(2),hDV) ;%  (APhi,Sv,mu,Pp,hDV)
            % find highest and lowest values of all 81 cases
            [minSHmaxVal,minSHmaxIdx]=min(AllSHMaxVals);
            [maxSHmaxVal,maxSHmaxIdx]=max(AllSHMaxVals);
            [minShminVal,minShminIdx]=min(AllShmaxVals);
            [maxShminVal,maxShminIdx]=max(AllShmaxVals);
            % find difference (delta) between middle value and
            % highest/lowest
            newSHmaxDetValue=mean([minSHmaxVal,maxSHmaxVal]);
            newShminDetValue=mean([minShminVal,maxShminVal]);
            deltaSHmax=abs([newSHmaxDetValue-minSHmaxVal]);% in PSI at the depth
            deltaShmin=abs([newShminDetValue-minShminVal]);
            vec(3)=deltaSHmax./dpth;
            vec(2)=deltaShmin./dpth;
            % change deterministic values slightly because deterministic
            % value isn't in middle of bounds necessarily
              disp(['changing SHmin value from ',num2str(hDV.data.stress.vals(2).*dpth),' to ',num2str(newShminDetValue),' a difference of ',num2str(diff([hDV.data.stress.vals(2).*dpth,newShminDetValue])),...
                  ' PSI because the deterministic value is not in the middle of the bounds'])
              disp(['changing SHmax value from ',num2str(hDV.data.stress.vals(3).*dpth),' to ',num2str(newSHmaxDetValue),' a difference of ',num2str(diff([hDV.data.stress.vals(3).*dpth,newSHmaxDetValue])),...
                  ' PSI because the deterministic value is not in the middle of the bounds'])                  
            % change deterministic SH and Sh
            hDV.data.stress.vals(2)=newShminDetValue./dpth;
            hDV.data.stress.vals(3)=newSHmaxDetValue./dpth;
%             
%             [SH,Sh] = getHorFromAPhi(APhi,Sv,mu,Pp,hDV)

        end
                   
        hDV.data.sigvals=vec; % store
        flag = checkProbabilisticEntries(hDV,noshow); % flag=1 if no issues
        if flag
            delete(figureDataEntryProb);
            modal(hfig,'on');
        end
    end
% check probabilistic data distribution bound entries
    function [flag] = checkProbabilisticEntries(hDV,noshow)
        flag=true; % flag=1 if no issues
        % by Rall
        % still need to check relative stress magnitudes etc, and if stress state
        % is possible.
        for kj=1:length(hDV.data.distributions.distriutionType);
            if hDV.data.distributions.vals(kj,1)<0
                edlgbox = errordlg(cat(2,'check ',hDV.data.distributions.distirbutionTxt{kj},' you entered a negative number: ',num2str(hDV.data.distributions.vals(kj,1)),...
                    'FSP will add and subtract a positive number from the deterministic value to get the bounds of the uniform distribution'));
                centerFigure(hDV.hfig,edlgbox);
                flag=false;
                
            end
        end
    end

    function convertAphiMuToStressVals(hDV)
        
        
    end

end




