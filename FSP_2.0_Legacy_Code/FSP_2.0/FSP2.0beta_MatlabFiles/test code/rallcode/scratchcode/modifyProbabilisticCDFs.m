
% change highlighting CDF curves for AGU talk
% make dashed grey lines CDF curves

% thisFault=hSV.plotdata.curfault(1);
thisFault=13;
indexes=1:hSV.data.fault.vals(1);
indexes(thisFault)=[];
set(hSV.plotdata.flinesprob(indexes),'color',[.5,.5,.5])
set(hSV.plotdata.flinesprob(:),'linestyle','-','linewidth',2)
set(hSV.plotdata.flinesprob(thisFault),'linewidth',3)

set(hSV.plotdata.flinesprob(indexes),'linestyle','--')






