function []=tornado()


% Names of the Y axis ticks
dat={'A-Map Event Rate' , .2190 , .2436 ;
       'b-value' , .18919 , .28915 ;
       'A-Map Spatial Distribution' , .20629 , .35622;
       'Maximum Magnitude' , .2053 , .2546;
       'GMPE sigma' , .18382 , .32121 ;
       'GMPE mean' , .2347 , .33403 ;
       'GMPE branch' , .1647 , .2792 ;
};

names = dat(:,1) ; 

% The base value is where the y axis is centered
base_value=0.2347;

low_vals  = cell2mat(dat(:,2)) ; 
high_vals = cell2mat(dat(:,3)) ;




% Sort the values based on the lower change
% Sort the higher values and the names arrays
%    using the same indices
[~ , ind]=sort(abs(low_vals-high_vals),'ascend');
high_vals=high_vals(ind);
low_vals=low_vals(ind);
names=names(ind);

% Create a figure and plot the low and high horizontally
figure
h = barh(high_vals,'r');
hold on
barh(low_vals)
bh = get(h,'BaseLine');
set(bh,'BaseValue',base_value);
title('Sensitivities','fontsize',14)
set(gca,'yticklabel',names)
set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
set(gca,'yticklabel',names)
set(gca,'fontsize',14)
xlabel('PGA [g]')
grid on

end
