function calldatmc_backedUp(src,evt,hDV)

title = 'Uniform Distribution bounds' ;
% i=6 ; % where does the data go, see dataentry.m
noshow = [4 9 11 12] ; %what not to show

if isfield(hDV.data,'distributions')
    vals = hDV.data.sigvals;
    if  isfield(hDV.data.stress.aphi,'sigvals') % sigma of aphi/mu not sigma of SH and Sh
        vals(2)=hDV.data.stress.aphi.sigvals(1);
        vals(3)=hDV.data.stress.aphi.sigvals(2);
    else
        hDV.data.stress.aphi.sigvals = [0,0];
    end
else
    vals = hDV.data.sigvals;  % backward compatible
    hDV.data.distributions.vals = vals ;
    
    % 0=constant and noshow, 1=uniform distribution, 2=constant, ...
    hDV.data.distributions.distriutionType=ones(12,1);
    hDV.data.distributions.distriutionType(noshow,1)=0;
end

thf=hDV.data.fault.thf  ;   %fault strikes
dips=hDV.data.fault.dipf ; %fault dips
SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
mufs=hDV.data.fault.muf ;
%     dpth=hDV.data.stress.vals(5);
sig = hDV.data.stress.vals(1:3)  ; sig=sort(sig,'descend');
ixSv = find(sig==hDV.data.stress.vals(1) ); %index of vertical stress
ixSv = ixSv(1) ; %if there are more than one
pp0 = hDV.data.stress.vals(6) ;
nu = hDV.data.reservoir.vals(5)  ;
biot = 1 ; %biot Coefficient
dp=0 ; %pressure perturbation


if all(thf==thf(1))
    strikeDetVal=thf(1);
else
    strikeDetVal= 'varying, ';% [cat(2,'between ',num2str(min(thf)),' and ',num2str(max(thf)))]
end
if all(dips==dips(1))
    dipsDetVal=dips(1);
else
    dipsDetVal= 'varying, ';% [cat(2,'between ',num2str(min(thf)),' and ',num2str(max(thf)))]
end
if all(mufs==mufs(1))
    mufsDetVal=mufs(1);
else
    mufsDetVal= ['varying, ',num2str(min(mufs)),' to ',num2str(max(mufs))];% [cat(2,'between ',num2str(min(thf)),' and ',num2str(max(thf)))]
end

% do we make this take aphi and mu? 
hDV.data.distributions.deterministicVals={hDV.data.stress.vals(1,1);hDV.data.stress.vals(2,1);hDV.data.stress.vals(3,1);ixSv;pp0;strikeDetVal;dipsDetVal;SHdir;dp;mufsDetVal;biot;nu};
% deterministicVals = hDV.data.distributions.deterministicVals;
deterministicValNums=nan(numel(hDV.data.distributions.deterministicVals),1);
for p=1:numel(hDV.data.distributions.deterministicVals)
    if ~ischar(hDV.data.distributions.deterministicVals{p,1})
deterministicValNums(p) = hDV.data.distributions.deterministicVals{p};
    end
end

txt2 = {'Vertical Stress Grad [',' psi/ft]' ;
    'Min Horiz. Grad [',' psi/ft]' ;
    'Max Horiz. Grad [',' psi/ft]' ;
    'ixSv: ' ,'';
    'Initial PP Grad [',' psi/ft]' ;
    'Strike Angles [',' degrees]' ;
    'Dip Angles [',' degrees]' ;
    'Max Horiz. Stress Dir [',' degrees]' ;
    'delta pore pressure','' ;
    'Friction Coeff Mu [',']' ;
    'Biot Parameter' ,'';
    'Poisson''s Ratio [' ,']'};

txt=cell(12,1);
for kj=1:12 % concatenate string with value from input
    %         txt2{kk,1},num2str(deterministicVals(kk)),txt2{kk,2}
    txt{kj,1}=cat(2,txt2{kj,1},num2str(hDV.data.distributions.deterministicVals{kj}),txt2{kj,2});
end

if hDV.data.stress.aphi.use >0; %not using APhi in GUI
    txt{2,1} = cat(2,'A Phi Parameter [',num2str(hDV.data.stress.aphi.vals(1)),']');
    txt{3,1} = cat(2,'reference mu [',num2str(hDV.data.stress.aphi.vals(2)),']');
    vals(2)=hDV.data.stress.aphi.sigvals(1);
    vals(3)=hDV.data.stress.aphi.sigvals(2);
    hDV.data.distributions.deterministicVals{2} = hDV.data.stress.aphi.vals(1);
    hDV.data.distributions.deterministicVals{3} = hDV.data.stress.aphi.vals(2);
    deterministicValNums(2)= hDV.data.stress.aphi.vals(1);
    deterministicValNums(3)=hDV.data.stress.aphi.vals(2);
end

% bound1OrPlusMinus2=[1;1;1;1;1;2;2;1;1;2;1;1]; % parameters that vary with each fault(2), or that vary with stresses(1)
% hDV.data.distributions.bound1OrPlusMinus2=bound1OrPlusMinus2;
dataentryProbabilistic(hDV,title,txt,vals,noshow,deterministicValNums);
hDV.data.distributions.distirbutionTxt=txt; % save text

%     indatacell ={sig',ixSv,pp0,thf,dips,SHdir,0,muf,biot,nu} ;
%     sigcell = indatacell ;
%     sigcell = {[.1 .1 .1],0,.05,5,5,5,0,.1,0,0}  ;

end


