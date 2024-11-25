
% debug Graphics:

% problem with some drivers: text of figures appears upside down
% toggle graphics renderer to try to debug


function toggleRenderer(src,~,hDV)
if nargin==0
    figurehande=gcf;
end
hDI.hfig=hDV.hfig;
set(hDV.plotdata.pint.hsf,'visible','on')
switch get(src,'value')
    case 1
        
        set(hDI.hfig,'Renderer','openGL');  % zbuffer, openGL, painters
        set(hDI.hfig,'RendererMode','manual');
        if ispc % if is on a pc, use software version of OpenGL renderer, not hardware, to avoid flipping text
            % had compatibility with some hardware with Doug - email 12/2016
            % -FRW
            opengl software
            disp('openGL software')
        else
            disp('not seen as a PC')
        end
    case 2
        
        set(hDI.hfig,'Renderer','openGL');  % zbuffer, openGL, painters
        set(hDI.hfig,'RendererMode','manual');
        if ispc % if is on a pc, use software version of OpenGL renderer, not hardware, to avoid flipping text
            % had compatibility with some hardware with Doug - email 12/2016
            % -FRW
            opengl hardware
            disp('openGL hardware')
        else
            disp('not seen as a PC')
        end
        
    case 3
        
        set(hDI.hfig,'Renderer','zbuffer');  % zbuffer, openGL, painters
        set(hDI.hfig,'RendererMode','manual');
        disp('renderer zbuffer - no opacity')
        set(hDV.plotdata.pint.hsf,'visible','off')
        
    case 4
        
        set(hDI.hfig,'RendererMode','manual');
        disp('renderermode  manual')
        
        
end























