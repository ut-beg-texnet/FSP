function [xs,ys] = setuparrowstressdirStereonet(hDV)

if nargin==0
    hDV = evalin('base','hSV');
end

% make SHmax arrow outside stereonet in Geomechanics tab
% by Rall Walsh
% make the stereonet axis a bit bigger in setupplotpanels:
% hDV.plotdata.pffot.ax3=axes('parent',pl,'position',[.55  .04  .45  .45],'fontsize',hDV.ftsz) ;

% this line can also go in setupplotpanels:
% set(findall(hDV.plotdata.pffot.ax3,'type','patch'),'visible','off')  % this prevents white patch objects from covering up stereonet poles if restacked


% point order for SHmax azimuth arrows
%               1
%               |
%               |
%           3   |   5
%            \  |  /
%             \ | /
%              2,4

sHAzimuth=hDV.data.stress.vals(4); %SH azimuth in degrees
trigAngle2=90-sHAzimuth; % angle form x axis
trigAngle=mod(trigAngle2+180,360); % add 180 degrees
shaftLength=0.3; 
bladeLength=0.15;
bladeHalfAngle=3; % degrees, from center of stereonet plus or minus to end of blade
radiusToArrowPoint=1.2;
axislims=[-1,1].* 1.5;

set(hDV.plotdata.pffot.ax3,'ylim',axislims,'xlim',axislims,'xtick',[],'ytick',[])
[XEdgeOfStereonet,YEdgeOfStereonet] = pol2cart(trigAngle*pi/180*[1 1],radiusToArrowPoint*[1,-1]); % arrow points
centerX=XEdgeOfStereonet(1); % center of points
centerY=YEdgeOfStereonet(1);
% calculate points
xs=[]; % clear
ys=[];
[xs(1),ys(1)] = pol2cart([trigAngle+bladeHalfAngle].*pi/180*[1],radiusToArrowPoint+bladeLength); % blade
xs(2)=centerX; % point
ys(2)=centerY;
[xs(3),ys(3)] = pol2cart([trigAngle]*pi/180*[1],radiusToArrowPoint+shaftLength); % shaft
xs(4)=centerX; % point
ys(4)=centerY;
[xs(5),ys(5)] = pol2cart([trigAngle-bladeHalfAngle]*pi/180*[1],radiusToArrowPoint+bladeLength); % blade
xs(6)=NaN; % disconnect 2 arrows
ys(6)=NaN;
centerX=XEdgeOfStereonet(2);  % second arrow point
centerY=YEdgeOfStereonet(2);
[xs(7),ys(7)] = pol2cart([mod(trigAngle+bladeHalfAngle,360)].*pi/180*[1],(radiusToArrowPoint+bladeLength)*-1); % blade: -1 reflects across center of stereonet
xs(8)=centerX;  % point
ys(8)=centerY;
[xs(9),ys(9)] = pol2cart([trigAngle]*pi/180*[1],(radiusToArrowPoint+shaftLength)*-1); % shaft
xs(10)=centerX;  %
ys(10)=centerY;
[xs(11),ys(11)] = pol2cart([mod(trigAngle-bladeHalfAngle,360)]*pi/180*[1],(radiusToArrowPoint+bladeLength)*-1); % blade

if nargin==0
    set(hDV.plotdata.sHmaxdir,'xdata',xs ,'ydata' , ys ,'color',[.5 .5 .5]) %max hor direction
end

end