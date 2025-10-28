function callbacksavesess(~,~,hDV)

datsviz=hDV.data;
% chkstring='DARREN PAIS URC 2016';  % changed in V0.99.1
chkstring='Saved FSP Session'; 
fname='FaultSlipPotential.mat';
uisave({'datsviz','chkstring'},fname)

end