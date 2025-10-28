classdef DASVdispQRA < hgsetget
    
    %--------------------------------------------------------------------------
    % (c) ExxonMobil Upstream Research Company, 2013++
    %     Drilling & Subsurface Technology Division
    %     Well Construction Section
    %  Authors: Pais, D; Payette, GS
    %
    % This is the display super-class. Comprises all basic properties and
    % methods inherited by all display windows.
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    %
    % Class properties
    %
    %--------------------------------------------------------------------------
    
    %----------------------------Public properties----------------------------%
    properties (Access = 'public')
        
        viz  = [] ; % figure visibility 'on' or 'off'
        hfig = [] ; % figure handle to main plotting figure
        name = [] ; % heading name for figure
        
    end
    
    %--------------------------------------------------------------------------
    %
    % Public methods
    %
    %--------------------------------------------------------------------------
    
    methods (Access = 'public')
        
        %----------------------------Class constructor----------------------------%
        function hDI = DASVdispQRA(name_in)
            hDI.viz ='on';
            hDI.name = name_in;
            [ll, bb, ww, hh] = hDI.setFigureSize();
            hDI.hfig = figure('tag',hDI.name,'visible',hDI.viz,'color',[1 1 1],...
                'units','pixels','position',[ll bb ww hh],'CloseRequestFcn',    ...
                {@exitDoublecheck,hDI},'tag','quit');
%                 @(src, event)hDI.delete()); 

            set(hDI.hfig,'MenuBar','none');
            set(hDI.hfig,'Name',hDI.name);
            set(hDI.hfig,'NumberTitle','off');
            set(hDI.hfig,'Renderer','openGL');  % zbuffer, openGL, painters
            set(hDI.hfig,'RendererMode','manual');
            
            
            set(hDI.hfig,'DefaultUicontrolUnits','normalized');
            set(hDI.hfig,'DefaultUicontrolFontsize',14);
            set(hDI.hfig,'PaperPositionMode','auto');
            if ispc % if is on a pc, use software version of OpenGL renderer, not hardware, to avoid flipping text
                % had compatibility with some hardware with Doug - email 12/2016
                % -FRW
                opengl software
            end
        end  % DASVdispQRA method (class constructor)
        
        
        %----------------------------Class destructor-----------------------------%
        function delete(hDI)
            
            if ishandle(hDI.hfig)
                delete(hDI.hfig);
            end
            
        end  % delete method (class destructor)
        
        %--------------------Constrain aspect ratio on figure---------------------%
        function mainfig_resize(hDI,pad)
            sz = get(hDI.hfig,'position');
            ss = get(0,'screensize');
            WS = ss(3);
            HS = ss(4);
            ll = pad * WS;
            bb = pad * HS;
            ww = (1 - 2*pad) * WS;
            hh = (1 - 2*pad) * HS;
            
            if ( (sz(3) < ww) || (sz(4) < hh) )
                sz = [ll bb ww hh];
            end
            
            set(hDI.hfig,'position',sz) ;
        end  % mainfig_resize method
        
        %--------------------------Set figure visibility--------------------------%
        function setviz(hDI,onoff)
            hDI.viz = onoff;
            set(hDI.hfig,'visible',hDI.viz) ;
        end  % setviz method
        
        %------------------------Toggle figure visibility-------------------------%
        function toggleviz(hDI)
            if strcmpi(hDI.viz,'on')
                setviz(hDI,'off');
            else
                setviz(hDI,'on');
            end
        end  % toggleviz method
        
        %---------------------------Set the figure name---------------------------%
        function setname(hDI,newname)
            
            hDI.name=newname;
            set(hDI.hfig,'Name',hDI.name);
            
        end  % setname method
        
    end  % PUBLIC methods
    
    %--------------------------------------------------------------------------
    %
    % Static methods
    %
    %--------------------------------------------------------------------------
    
    methods (Static)
        
        %-----------------------------Set figure size-----------------------------%
        function [ll, bb, ww, hh] = setFigureSize()
            
            Apref=0.8*1920; Bpref=0.8*1200; ARpref=Apref/Bpref; %prefered size at build
            ss  = get(0,'screensize'); A=ss(3); B=ss(4);
            
            if A>Apref && B>Bpref
                out = [(A-Apref)/2 (B-Bpref)/2 Apref Bpref] ;
            elseif A/B>=ARpref
                out = [(A-ARpref*B)/2 0 ARpref*B B] ;
            elseif  A/B<ARpref
                out = [0 (B-A/ARpref)/2 A A/ARpref] ;
            else
                error('No match for screen resoulution') ;
            end
            
            %    ss  = get(0,'screensize');
            %    WS  = ss(3);
            %    HS  = ss(4);
            %    pad = 0.2; %scales figure by this amount and centers it
            %    ll  = pad*WS;
            %    bb  = pad*HS;
            %    ww  = (1-2*pad)*WS;
            %    hh  = (1-2*pad)*HS ;
            
            ll=out(1) ; bb=out(2) ; ww=out(3) ; hh=out(4) ;
            
        end  % setFigureSize method
        
        %----------------------------Display messages-----------------------------%
        function disp_msg(msg,outtype)
            switch outtype
                case 'screen'
                    disp([datestr(now),' - ',msg]) ;
                    
            end
        end  % disp_msg method
        
    end  % STATIC methods
    
end % classdef

%--------------------------------------------------------------------------
%
% End of class: class dismissed
%
%--------------------------------------------------------------------------