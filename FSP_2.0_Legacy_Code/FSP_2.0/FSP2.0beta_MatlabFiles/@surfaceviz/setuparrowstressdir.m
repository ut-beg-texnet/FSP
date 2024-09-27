function setuparrowstressdir(hDV)

% point order for SHmax azimuth arrows
%               1
%               |
%               |
%           5   |   7
%            \  |  /
%             \ | /
%              3,6
%             / |\
%            /  | \
%           8   |  4
%               |
%               |
%               2

sHAzimuth=hDV.data.stress.vals(4); %SH azimuth in degrees
trigAngle=90-sHAzimuth; % angle form x axis
trigAngle2=mod(trigAngle+180,360); % add 180 degrees

shaftLength=0.9;
bladeLength=0.3;
bladeHalfAngle=14; % degrees
centerX=0; % center of points
centerY=0;
% calculate points
xs=[];
ys=[];
ys(1)=centerX+shaftLength.*sind(trigAngle);
xs(1)=centerX+shaftLength.*cosd(trigAngle);
ys(2)=centerY+shaftLength.*sind(trigAngle2);
xs(2)=centerY+shaftLength.*cosd(trigAngle2);
xs(3)=centerX;
ys(3)=centerY;
xs(4)=centerX+bladeLength.*cosd(trigAngle2+bladeHalfAngle);
ys(4)=centerY+bladeLength.*sind(trigAngle2+bladeHalfAngle);
xs(5)=centerX+bladeLength.*cosd(trigAngle+bladeHalfAngle);
ys(5)=centerY+bladeLength.*sind(trigAngle+bladeHalfAngle);
xs(6)=centerX;
ys(6)=centerY;
xs(7)=centerX+bladeLength.*cosd(trigAngle-bladeHalfAngle);
ys(7)=centerY+bladeLength.*sind(trigAngle-bladeHalfAngle);
xs(8)=centerX+bladeLength.*cosd(trigAngle2-bladeHalfAngle);
ys(8)=centerY+bladeLength.*sind(trigAngle2-bladeHalfAngle);

set(hDV.plotdata.inputMap.hPlotArrowSH,'xdata',xs,'ydata',ys)
set(hDV.plotdata.inputMap.ax3,'xlim',[-1,1],'ylim',[-1,1],'xtick',[],'ytick',[]) ;

%fix aspect ratio to match the aspect ratio of faults plot
daspect(hDV.plotdata.inputMap.ax3,pbaspect(hDV.plotdata.inputMap.ax1));


end