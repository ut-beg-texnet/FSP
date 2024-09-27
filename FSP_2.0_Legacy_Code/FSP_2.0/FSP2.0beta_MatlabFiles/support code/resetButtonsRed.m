
% By Rall
% Stanford
% after new data loaded in, turn buttons red, so that they need to be
% calculated. 

% if just hydrologic data changed, no need to turn geomechanics buttons
% red. 

% therefore buttonIndexToReset specifies which buttons to reset. 

% just inputting hDV resets all buttons that can be red. 

% button index is workflow order: ie calculate first, prob geomech
% second,prob hyd third, summary 4th. 

% Replaces these lines
%                     set(hDV.plotdata.pprob.hbutMC,'backgroundcolor',hDV.colors.red)
%                     set(hDV.bCalc,'backgroundcolor',hDV.colors.red)
%                     set(hDV.plotdata.hydprob.hbutMC,'backgroundcolor',hDV.colors.red)
%                     set(hDV.plotdata.pprob.summaryPlotButton,'backgroundcolor',hDV.colors.red)

function  resetButtonsRed(hDV,buttonIndexToReset)


                 % calculate button, Monte Carlo Geomech     Prob. hydrology button       Sumary plot button
buttonsInWorkflowOrder= [hDV.bCalc;hDV.plotdata.pprob.hbutMC;hDV.plotdata.hydprob.hbutMC;hDV.plotdata.pprob.summaryPlotButton];
if nargin ==1 % if no first button specified, turn all red
    buttonIndexToReset=1:length(buttonsInWorkflowOrder);
end

set(buttonsInWorkflowOrder(buttonIndexToReset),'backgroundcolor',hDV.colors.red)

end