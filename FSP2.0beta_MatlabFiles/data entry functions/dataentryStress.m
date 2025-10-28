% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface
%
% Modified by Suvrat Lele   March 2018  for FSP 2.0
% Subsurface Mechanics and Induced Seismicity, D&S

%GUI for entry of in situ stress data

function valsout = dataentryStress(hDV,title,txt,vals,i)

hfig = hDV.hfig;
modal(hfig,'off');

hpos = get(hfig,'position') ;
figureStressEntry = figure('tag',title,'color',[1 1 1],'units','pixels','visible','off');
pos = [0 0 .6*hpos(3) .7*hpos(4)];
set(figureStressEntry,'position',pos) ;
set(figureStressEntry,'MenuBar','none');
set(figureStressEntry,'Name',title);
set(figureStressEntry,'NumberTitle','off');
set(figureStressEntry,'DefaultUicontrolUnits','normalized');
set(figureStressEntry,'DefaultUicontrolFontsize',14);
set(figureStressEntry,'PaperPositionMode','auto');
set(figureStressEntry,'closerequestfcn',@crf) ;
centerFigure(hfig,figureStressEntry);
set(figureStressEntry,'visible','on') ; 



% Show the radio toggle for horizontal stress model, gradients or A-Phi
%
radio1 = uicontrol('Style','radio','String','Specify All Three Stress Gradients [psi/ft]', 'Tag', 'strModel_gradients', ...
    'pos',[.1 .88 .4 .04],'parent',figureStressEntry,'BackgroundColor',[1 1 1],'HandleVisibility','on', 'value',0,...
    'ToolTipString','Directly specifying vertical and the two horizontal stress gradients, if available, is the best FSP stress model option');
radio2 = uicontrol('Style','radio','String','Use A-Phi Model','Tag', 'strModel_aphi', ...
    'pos',[.1 .8 .4 .04],'parent',figureStressEntry,'BackgroundColor',[1 1 1],'HandleVisibility','on', 'value', 0,...
    'ToolTipString','If horizontal stress data are not available then these may be estimated using the A-Phi model');


% Create text boxes for entering stress related values + 2 for A-Phi and 
% Posson's ratio
%
% Initially set visible off until a stress model is selected by user using
% above radio buttons


pos = [0.05 0 .45 .04] ;  c2=[.5 0 -.15 0] ;
valtxtH = zeros(length(vals)+2,1) ;
txtH = zeros(length(vals)+2,1) ;

% Vertical stress
%
txtH(1) = uicontrol('parent',figureStressEntry,'style','text','string',txt{1},...
    'pos',pos+[0 .70 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(1)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(1)),...
     'pos',c2+pos+[0 .70 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');

% Horizontal MAX stress
%
txtH(3) = uicontrol('parent',figureStressEntry,'style','text','string',txt{3},...
    'pos',pos+[0 .63 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(3)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(3)),...
     'pos',c2+pos+[0 .63 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');

% Horizontal MIN stress
%
txtH(2) = uicontrol('parent',figureStressEntry,'style','text','string',txt{2},...
    'pos',pos+[0 .56 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(2)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(2)),...
     'pos',c2+pos+[0 .56 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');
 
% Horizontal max stress azimuth
%
txtH(4) = uicontrol('parent',figureStressEntry,'style','text','string',txt{4},...
    'pos',pos+[0 .46 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(4)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(4)),...
     'pos',c2+pos+[0 .46 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');
 
% Pore pressure gradient
%
txtH(6) = uicontrol('parent',figureStressEntry,'style','text','string',txt{6},...
    'pos',pos+[0 .39 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(6)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(6)),...
     'pos',c2+pos+[0 .39 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');
 
% Reference depth
%
txtH(5) = uicontrol('parent',figureStressEntry,'style','text','string',txt{5},...
    'pos',pos+[0 .32 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off');     
valtxtH(5)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(vals(5)),...
     'pos',c2+pos+[0 .32 0 0],'fontsize',12,'fontunits','normalized',...
     'backgroundcolor',[1 1 1], 'visible', 'off');


 
PoissonVisible = 'off';   % SET ON FOR INTERNAL EXXONMOBIL - XTO VERSION ONLY <<<<<-----------

 
% Add input for Poisson's ratio on stress data entry tab; it is stored in
% reservoir dataset -- not changed for compatibility with saved sessions
% from older FSP versions
%
poissonTip = sprintf('Use default value of 0.5 for analyzing FSP in basement or any layer other than the injection layer\nUse a value appropriate for injection layer rock for analysing FSP in that layer');
txtH(8) = uicontrol('parent',figureStressEntry,'style','text','string',hDV.data.reservoir.txt(5),...
    'pos',pos+[0 .22 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off', 'ToolTipString', poissonTip);
%
valtxtH(8)= uicontrol('parent',figureStressEntry,'style','edit','string',num2str(hDV.data.reservoir.vals(5)),...
    'pos',c2+pos+[0 .22 0 0],'fontsize',12,'fontunits','normalized',...
    'backgroundcolor',[1 1 1], 'visible', 'off', 'ToolTipString', poissonTip);



%these are for the APhi and mu inputs
aphivals=hDV.data.stress.aphi.vals; 
txtH(7) = uicontrol('parent',figureStressEntry,'style','text','string','A-Phi Parameter',...
    'pos',pos+[0 .63 0 0],'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1], 'visible', 'off', ...
    'ToolTipString','A-Phi values for various regions may be found in published papers, e.g. Lund-Snee and Zoback (2018) for the Permian Basin');
valtxtH(7) = uicontrol('parent',figureStressEntry,'style','edit','string',num2str(aphivals(1)),...
    'pos',c2+pos+[0 .63 0 0],'fontsize',12,'fontunits','normalized','backgroundcolor',[1 1 1],'visible','off',...
    'ToolTipString','A-Phi values for various regions may be found in published papers, e.g. Lund-Snee and Zoback (2018) for the Permian Basin');


checkbox_modAphi = uicontrol('parent',figureStressEntry,'Style','checkbox',...
    'String','Min Horiz Stress Grad Available [psi/ft]', 'Tag', 'modified_Aphi', ...
    'pos',pos+[0.16 .56 -0.13 0],'fontsize',12,'BackgroundColor',[1 1 1],...
    'HandleVisibility','on', 'value',0,'visible','off',...
    'ToolTipString', 'Providing minimum horizontal stress value, if available, increases accuracy of the stress model');

% If using A Phi model based on prior input, set modified A-Phi checkbox
% and shmin box in correct state
%
if hDV.data.stress.aphi.use == 12
    set(checkbox_modAphi,'value',1);  % set Modified A-Phi checkbox to checked state
    if get(radio2,'value')
        set(valtxtH(2),'visible','on'); % make text box for shmin visible if A-Phi model selected by user
    end    
end



% Set callback function for the above stress model buttons 
%
set(radio1,'Callback',{@setStrModel,radio1,radio2,txtH,valtxtH,checkbox_modAphi});
set(radio2,'Callback',{@setStrModel,radio1,radio2,txtH,valtxtH,checkbox_modAphi});

% Set callback funciton for modified A-Phi checkbox
%
set(checkbox_modAphi,'Callback',{@setModAphi,radio2,valtxtH});



%
% Configure all input boxes based on previously entered/loaded input stress 
% data, if any
%

% If aphi.use in the data structure is nonzero then data for A-phi model
% has been specified earlier
%
if hDV.data.stress.aphi.use > 0
    
    % Set radio button for A-Phi model and force callback (as if the radio
    % button was clicked by user)
    set(radio2,'value',1);
    setStrModel(radio2,[],radio1,radio2,txtH,valtxtH,checkbox_modAphi);

% Otherwise, if aphi.use==0 AND sHmax >0, then gradients model data has been
% specified earlier; since all values are initizied to zero at beginning, 
% just checking aphi.use==0 won't work
%
elseif hDV.data.stress.aphi.use == 0 && hDV.data.stress.vals(3) > 0
    
    % Set radio button for gradients model and force callback (as if the radio
    % button was clicked by user)
    set(radio1,'value',1);
    setStrModel(radio1,[],radio1,radio2,txtH,valtxtH,checkbox_modAphi);
    
end



% Create OK button at the bottom of the screen
%
uicontrol('parent',figureStressEntry,'style','pushbutton','string','OK',...
    'pos',[.32 .05 .36 .07],'fontsize',12,'fontunits','normalized',...
    'callback',{@butcall,hDV,valtxtH,radio1,radio2,checkbox_modAphi,i,figureStressEntry});


% Return value for dataentryStress 
%
valsout=0;



    % Callback function for stress model radio buttons
    %
    function setStrModel(src,~,radio1,radio2,txtH,valtxtH,checkbox_modAphi)
        
        switch get(src,'Tag')
            case 'strModel_gradients'
                if get(radio1,'value')
                    % User selected strModel_gradients here
                    % Deselect strModel_aphi
                    set(radio2,'value',0);
                else
                    % User deselected strModel_gradients 
                    % strModel_aphi should already be deselected; i.e. no model selected
                end
            case 'strModel_aphi'
                if get(radio2,'value')
                    % User selected strModel_aphi here
                    % Deselect strModel_gradients
                    set(radio1,'value',0);
                else
                    % User deselected strModel_aphi 
                    % strModel_gradients should already be deselected; i.e. no model selected
                end
        end
        
        % If gradient model is selected then show appropriate input boxes
        %
        if get(radio1,'value')
            set(txtH(1),'visible','on');
            set(txtH(2),'visible','on');
            set(txtH(3),'visible','on');
            set(txtH(4),'visible','on');
            set(txtH(5),'visible','on');
            set(txtH(6),'visible','on');
            set(txtH(7),'visible','off');
            set(txtH(8),'visible',PoissonVisible);
            set(checkbox_modAphi,'visible','off');
            set(valtxtH(1),'visible','on');
            set(valtxtH(2),'visible','on');
            set(valtxtH(3),'visible','on');
            set(valtxtH(4),'visible','on');
            set(valtxtH(5),'visible','on');
            set(valtxtH(6),'visible','on');
            set(valtxtH(7),'visible','off');
            set(valtxtH(8),'visible',PoissonVisible);
        end
        
        if get(radio2,'value')
            set(txtH(1),'visible','on');
            set(txtH(2),'visible','off');
            set(txtH(3),'visible','off');
            set(txtH(4),'visible','on');
            set(txtH(5),'visible','on');
            set(txtH(6),'visible','on');
            set(txtH(7),'visible','on');
            set(txtH(8),'visible',PoissonVisible);
            set(checkbox_modAphi,'visible','on');
            set(valtxtH(1),'visible','on');
            if get(checkbox_modAphi,'value')
                set(valtxtH(2),'visible','on');
            else
                set(valtxtH(2),'visible','off');
            end
            set(valtxtH(3),'visible','off');
            set(valtxtH(4),'visible','on');
            set(valtxtH(5),'visible','on');
            set(valtxtH(6),'visible','on');
            set(valtxtH(7),'visible','on');
            set(valtxtH(8),'visible',PoissonVisible);
        end
        
        if ~(get(radio1,'value') || get(radio2,'value'))
            set(txtH(1),'visible','off');
            set(txtH(2),'visible','off');
            set(txtH(3),'visible','off');
            set(txtH(4),'visible','off');
            set(txtH(5),'visible','off');
            set(txtH(6),'visible','off');
            set(txtH(7),'visible','off');
            set(txtH(8),'visible','off');
            set(checkbox_modAphi,'visible','off');
            set(valtxtH(1),'visible','off');
            set(valtxtH(2),'visible','off');
            set(valtxtH(3),'visible','off');
            set(valtxtH(4),'visible','off');
            set(valtxtH(5),'visible','off');
            set(valtxtH(6),'visible','off');
            set(valtxtH(7),'visible','off');
            set(valtxtH(8),'visible','off');
        end
        
    end % callback function for radio buttons


    % Callback function for Modified A-Phi model checkbox
    %
    function setModAphi(src,~,radio2,valtxtH)
        
        if get(radio2,'value')  % A-Phi radio button is selected
            if get(src,'value')  
                % this function called because checkbox for Modified A-Phi is checked by user
                set(valtxtH(2),'visible','on'); % make text box for shmin visible if A-Phi model selected by user
            else
                % this function called because checkbox for Modified A-Phi is unchecked by user
                set(valtxtH(2),'visible','off');
            end
        end
        
    end % callback functino for checkbox


    function crf(~,~)
        hDV.tempvec=[];
        delete(gcf);
        modal(hfig,'on');
    end


    %this must tie back to callbackdatabuts.m (SPL: does this still apply?)
    % Callback function for OK button
    %
    function butcall(~,~,hDV,htxt,radio1,radio2,checkbox_modAphi,i,figureStressEntry)
        
        vec=zeros(length(htxt)-2,1);  % last two, A-phi and nu not part of stress data structure
        for n=1:1:length(vec)
            vec(n)=str2double(get(htxt(n),'string'));
        end
        
        APhi = str2double(get(htxt(7),'string')); % A-Phi parameter input
        nu = str2double(get(htxt(8),'string')); % Poisson's ratio input
        
        % Separate stress variables, created to improve readability of code below 
        %
        Sv = vec(1);
        Sh = vec(2);
        SH = vec(3);
        
        aphi_use = [];
        
        if get(radio1,'Value')
            disp('Using input horizontal stresses') ;
            aphi_use = 0; 
            
            % Check that stress magnitudes and order ok
            %
            if Sv<=0 || Sh<=0 || SH<=0
                edlgBox=errordlg('Stress gradient values should be > 0','Data Error') ; centerFigure(hDV.hfig,edlgBox);
            elseif Sh>SH
                edlgBox=errordlg('Minimum Horiz Stress should be lower than Maximum Horiz Stress','Data Error') ; centerFigure(hDV.hfig,edlgBox);
            end

        elseif get(radio2,'Value')
                      
            if ~get(checkbox_modAphi,'Value')
                disp('Using APhi for horiz stresses; NOW WITH REF_MU EQUAL TO FAULT FRIC COEFF') ;
                aphi_use = 11;
            else 
                disp('Using MODIFIED APhi model') ;
                aphi_use = 12;
            end
            
            if ~(APhi>=0 && APhi<=3)
                edlgBox= errordlg(['Check Aphi [0,3] and mu (>0) range. You have mu=',num2str(mu),' and APhi =',num2str(APhi)]) ; centerFigure(hDV.hfig,edlgBox);
            end
            
        else
            % No stress model selected
            edlgBox= errordlg('Please select a stress model') ; centerFigure(hDV.hfig,edlgBox);
        end
        
        
        % Set buttons red if any values were changed
        %
        if ~(all(hDV.data.stress.vals==vec) && hDV.data.reservoir.vals(5)==nu && hDV.data.stress.aphi.use==aphi_use && hDV.data.stress.aphi.vals(1) == APhi)
            % Set Calculate and geomechanics button red 
            resetButtonsRed(hDV) % turn all buttons red
        end
        
        % Calculate A-Phi value based on gradients, or stresses based on
        % A-Phi model so that those boxes are pre-populated next time data
        % input window is opened (similar to behavior in FSP 1.0)
        %
        if aphi_use==0  
            % Gradient model is being used; calcualte corresponding Aphi
            % Note that only APhi is calculated since ref_mu is no longer used
            % 
            if Sv > SH  % normal faulting
                APhi = (SH-Sh)/(Sv-Sh);
            elseif Sv<Sh  % Reverse faulting
                APhi = 2 + (Sh-Sv)/(SH-Sv);
            else
                APhi = 2 - (Sv-Sh)/(SH-Sh);
            end
        else
            % APhi model used; calculate corresponding stresses
            %
            [SH,Sh] = getHorFromAPhi(APhi,hDV.data.fault.vals(2),Sh,Sv,vec(6),aphi_use);
            vec(2) = Sh;
            vec(3) = SH;
        end
        
        
        % Store all input data in the data structure
        %
        hDV.data.stress.vals=vec;
        hDV.data.reservoir.vals(5)=nu; % Update Poisson's ratio in reservoir dataset
        hDV.data.stress.aphi.use     = aphi_use;
        hDV.data.stress.aphi.vals(1) = APhi;
        hDV.data.stress.aphi.vals(2) = NaN; % No separate ref mu input now
        
        
        delete(figureStressEntry);
        modal(hfig,'on');
        
    end % callback function for OK button


end




