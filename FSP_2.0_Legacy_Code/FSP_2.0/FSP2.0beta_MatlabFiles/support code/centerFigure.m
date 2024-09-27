function centerFigure(hfig1,hfig2)

%--------------------------------------------------------------------------
% (c) ExxonMobil Upstream Research Company, 2013++ 
%     Drilling & Subsurface Technology Division
%     Well Construction Section
%  Authors: Payette, GS
%
%  This function centers the position of hfig2 based on the position of
%  hfig1
%--------------------------------------------------------------------------
set(hfig2,'units',get(hfig1,'units')) % make sure figures have same units - FRW
pos1    = get(hfig1,'position');
pos2    = get(hfig2,'position');
pos2(1) = pos1(1) + 0.5 * (pos1(3) - pos2(3));
pos2(2) = pos1(2) + 0.5 * (pos1(4) - pos2(4));
set(hfig2,'visible','on','position',pos2)

%--------------------------------------------------------------------------
% 
% End of function
% 
%--------------------------------------------------------------------------
