% ExxonMobil Upstream Research Company 2016
% Darren Pais
% Surveillance and Optimization, Reservoir 
% Induced Seisimicty Integrated Team, Drilling and Subsurface

%Callback for running monte carlo Geomechanics

function callrunmc(src,~,hDV)

nfaults = hDV.data.fault.vals(1) ;

%get sensitivity bar values (needs to match setupplotpanels!!!)
hDV.plotdata.results.barlowdata=cell(nfaults,1);
hDV.plotdata.results.barhighdata=cell(nfaults,1);
hDV.plotdata.results.barnom=cell(nfaults,1);

%wait bar
h = waitbar(0,'Calculating ...','Name','Monte-Carlo Engine');
centerFigure(hDV.hfig,h)


outs=cell(nfaults,1); 
     maxData=1; % psi Xlimit   

for k=1:1:nfaults
    %check for close waitbar
    if ~ishandle(h)
        disp('Monte Carlo Premature Shutdown') ; 
        break
    end

    thf=hDV.data.fault.thf(k) ;   %fault strikes
    dips=hDV.data.fault.dipf(k); %fault dips
    SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
    muf=hDV.data.fault.muf(k);
    dpth=hDV.data.stress.vals(5);
    sig = hDV.data.stress.vals(1:3)*dpth;

    pp0 = hDV.data.stress.vals(6)*dpth ;
    nu = hDV.data.reservoir.vals(5)  ;
    biot = 1 ; %biot Coefficient
    dp=0 ; %pressure perturbation
    indatacell ={sig',0.00,pp0,thf,dips,SHdir,dp,muf,biot,nu} ;
    
    %default to zero variability
    if ~isfield(hDV.data,'sigvals')
        hDV.data.sigvals=zeros(12,1);
    end
    v=hDV.data.sigvals';
    
    v_sig = v(1:3); % stress sig
    sigcell = {v_sig*dpth,v(4),v(5)*dpth,...
        v(6),v(7),v(8),v(9),v(10),v(11),v(12)}  ; %careful about stress units
    
    if hDV.data.stress.aphi.use
        % Input 11th cell entry for aphi
        indatacell{11} = hDV.data.stress.aphi.vals(1);
        if ~isfield(hDV.data.stress.aphi,'sigvals');% legacy
            hDV.data.stress.aphi.sigvals = [0,0];
        end        
        sigcell{11} = hDV.data.stress.aphi.sigvals(1);    
    end

    
    % MONTE CARLO engine run
    %
    if isfield(hDV.data,'nrunsGeomech')
        nrunsGeomech=hDV.data.nrunsGeomech;
    else
        nrunsGeomech=1000;
    end
    %
    [outs{k},hDV.data.distributionsinData{k,1}] = monte_carlo(@mohrs_3D,indatacell,sigcell,hDV,nrunsGeomech) ;
    [f,x] = ecdf(outs{k}) ;
    
    % Display negative pressures to slip as zero
    %
    x(x<0) = 0; 
    
    % Update delpa_pp probability lines data in data structure (plots will update automatically)
    %
    set(hDV.plotdata.flinesprob(k),'Xdata',x,'Ydata',f) ;
    
    
    
    %capture the distribution bounds for of the data
    lev0=0.01;lev1=.2; lev2=.5 ; 
    hDV.plotdata.results.nomlow{k}= max(x(f<lev1));
    hDV.plotdata.results.nomhigh{k}= max(x(f<lev2));
    hDV.plotdata.results.nomBottom{k}= max(x(f<lev0));
    
    waitbar(k / nfaults,h,['Calculating ',num2str(100*k / nfaults,'%3.0f'), '% ...']) ;
    
    
    ptile=0.16; %corresponding to +/- 1 sigma when using normal distribution .. not used for uniform distrubution
    [outup , outdown , nom] = sensitivity(@mohrs_3D,indatacell,sigcell,hDV,[]);
    
    % Store in data structure for plotting
    %
    hDV.plotdata.results.barnom{k}=nom;       
       
    if hDV.data.stress.aphi.use==11 % then input 11th entry, aphi
        hDV.plotdata.results.barlowdata{k}=[outdown{1}(1) outdown{3} outdown{4} outdown{5} outdown{6} outdown{8} outdown{11}];
        hDV.plotdata.results.barhighdata{k}=[outup{1}(1) outup{3} outup{4} outup{5} outup{6} outup{8}   outup{11}];
        
    elseif hDV.data.stress.aphi.use==12 % additionally add shmin
        hDV.plotdata.results.barlowdata{k}=[outdown{1}(1:2) outdown{3} outdown{4} outdown{5} outdown{6} outdown{8} outdown{11}];
        hDV.plotdata.results.barhighdata{k}=[outup{1}(1:2) outup{3} outup{4} outup{5} outup{6} outup{8}   outup{11}];
        
    elseif hDV.data.stress.aphi.use==0
        hDV.plotdata.results.barlowdata{k}=[outdown{1} outdown{3} outdown{4} outdown{5} outdown{6} outdown{8}];
        hDV.plotdata.results.barhighdata{k}=[outup{1} outup{3} outup{4} outup{5} outup{6} outup{8}];
    
    end
    
    % Display negative pressures to slip as zero - update above stored data
    %
    if hDV.plotdata.results.barnom{k}<0; hDV.plotdata.results.barnom{k}=0; end
    hDV.plotdata.results.barlowdata{k}(hDV.plotdata.results.barlowdata{k}<0) = 0;
    hDV.plotdata.results.barhighdata{k}(hDV.plotdata.results.barhighdata{k}<0) = 0;
    
    
    % find maximum Xlimit value
    maxData=max([maxData,max(x)]);
    
end % loop over faults



% variable names for tornado chart
if hDV.data.stress.aphi.use ==11
        hDV.plotdata.pprob.names = {'Vert Stress Grad' ;
                 'Pore Press Grad ';
                 'Strike of Fault ' ;
                 'Dip of Fault    ' ;
                 'SHmax Azimuth   ' ;
                 'Friction Coeff  ';
                 'A-Phi           ';};
elseif hDV.data.stress.aphi.use ==12
    hDV.plotdata.pprob.names = {'Vert Stress Grad' ;
                 'Shmin Gradient  ';
                 'Pore Press Grad ';
                 'Strike of Fault ' ;
                 'Dip of Fault    ' ;
                 'SHmax Azimuth   ' ;
                 'Friction Coeff  ';
                 'A-Phi           ';};
else
      hDV.plotdata.pprob.names = {'Vert Stress Grad' ;
                 'Shmin Gradient  ';
                 'SHmax Gradient  ';
                 'Pore Press Grad ';
                 'Strike of Fault ' ;
                 'Dip of Fault    ' ;
                 'SHmax Azimuth   ' ;
                 'Friction Coeff  '};
end



% if geomechanics tab hasn't been run
if ~ isfield(hDV.plotdata.results,'outs')
    calcengine(hDV,'GEOMECHANICS')
end

%set the colors on plot
setprobcolor(hDV);

%refresh colorbar
%callcbarPtxt([],[],hDV);

       
%empty remaining
set(hDV.plotdata.flinesprob(nfaults+1:end),'xdata', NaN ,'ydata' , NaN) ;
%set(hDV.plotdata.flinesprobtxt(nfaults+1:end),'position', [NaN NaN])

%xlim(hDV.plotdata.paly.ax2,[0 maxdp]) ;
ylim(hDV.plotdata.pprob.ax2,[0 1]) ;

%remove wait bar and change color of button to green
if ishandle(h)
    delete(h) ;
    set(src,'backgroundcolor',hDV.colors.green); 
end 

%refresh panel
refreshplotdata(hDV,hDV.currtab.name); %plot data for selected panel refreshed
 

end


