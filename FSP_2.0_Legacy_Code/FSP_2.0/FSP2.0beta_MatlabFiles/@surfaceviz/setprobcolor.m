function setprobcolor(hDV)
cv=hDV.plotdata.results.outs.ppfail; nfaults = hDV.data.fault.vals(1) ; 
for j=1:1:nfaults
    cl = getcolor(flipud(hDV.cmapGYR),cv(j),hDV.plotdata.mincfc,hDV.plotdata.maxcfc) ;
    set(hDV.plotdata.flinesprob(j),'color',cl);
end

end
