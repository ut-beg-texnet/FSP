

% By Rall Walsh
% Stanford, August 2016
% manually add tick lines to stereonet polar plot at various dip angles
% still need to check equal area vs equal angle

% add these lines in setupplotpanels to run: 
%              [X5,Y5]=addStereonetLines(hDV);
%             plot(hDV.plotdata.pffot.ax3,X5,Y5,'k:'); % make sure this is correct - equal area vs equal angle? 

function [X,Y]=addStereonetLines(hDV)

if nargin==0
    
    figure
    polar(0,1)
    
end

RStereonet=1;
% radii to plot
dipLinesDegree=15:15:75;
radii=cosd(dipLinesDegree).*RStereonet;
theta=0:pi/180:2*pi;
[allRs,allThetas]=meshgrid(radii,theta);
allRs=[allRs;NaN(1,length(dipLinesDegree))];
allThetas=[allThetas;NaN(1,length(dipLinesDegree))];
 [X,Y]=pol2cart(allThetas,allRs);

if nargin == 0

    plot(X,Y,'k:')
end
end