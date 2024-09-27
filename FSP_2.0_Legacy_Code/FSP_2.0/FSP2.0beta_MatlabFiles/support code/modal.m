function modal(hfig,onoff)

%--------------------------------------------------------------------------
% (c) ExxonMobil Upstream Research Company, 2014++ 
%     Drilling & Subsurface Technology Division
%     Well Construction Section
%  Authors: Pais, D; Payette, GS
%
% Makes a figure modal or turns modal off
%--------------------------------------------------------------------------

if ( strcmpi(get(hfig,'visible'),'on') )
   
   sd.jFigPeer = get(handle(hfig),'JavaFrame');
   sd.jWindow  = sd.jFigPeer.fHG1Client.getWindow;
   
   switch onoff
      case 'on'
         sd.jWindow.setEnabled(true);
      case 'off'
         sd.jWindow.setEnabled(false);
   end
   
end

%--------------------------------------------------------------------------
% 
% End of file
% 
%--------------------------------------------------------------------------