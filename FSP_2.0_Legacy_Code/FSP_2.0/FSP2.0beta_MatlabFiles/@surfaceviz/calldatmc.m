%
% Modified by Suvrat Lele for FSP 2.0 March 2018; changes related to A-Phi
% Subsurface Mechanics and Induced Seismicity, Drilling & Subsurface
% ExxonMobil URC

function calldatmc(src,evt,hDV)

title = 'Uniform Distribution bounds' ;
% i=6 ; % where does the data go, see dataentry.m
noshow = [4 9 11 12] ; %what not to show
vals = hDV.data.sigvals;

% Set variability for noshow variables to zero (so that any values not
% visible to user do not change anything
hDV.data.sigvals(4) = 0;
hDV.data.sigvals(9) = 0;
hDV.data.sigvals(11) = 0;
hDV.data.sigvals(12) = 0;
        

if  isfield(hDV.data.stress.aphi,'sigvals') % sigma of aphi/mu not sigma of SH and Sh
    if hDV.data.stress.aphi.use
        vals(13)=hDV.data.stress.aphi.sigvals(1);
    end
else
    hDV.data.stress.aphi.sigvals = [0,0];
end

if ~isfield(hDV.data,'distributions')
    hDV.data.distributions.vals = vals ;
    % 0=constant and noshow, 1=uniform distribution, 2=constant, ...
    hDV.data.distributions.distriutionType=ones(length(vals),1);
    hDV.data.distributions.distriutionType(noshow,1)=0;
end

thf=hDV.data.fault.thf  ;   %fault strikes
dips=hDV.data.fault.dipf ; %fault dips
SHdir = hDV.data.stress.vals(4) ; %max horiz stress direction
mufs=hDV.data.fault.muf ;
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
hDV.data.distributions.deterministicVals={hDV.data.stress.vals(1,1);hDV.data.stress.vals(2,1);hDV.data.stress.vals(3,1);0.00;pp0;strikeDetVal;dipsDetVal;SHdir;dp;mufsDetVal;biot;nu};
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
    'IGNORE PARAMETER4: ' ,'';
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

if hDV.data.stress.aphi.use
    txt{13,1} = cat(2,'A Phi Parameter [',num2str(hDV.data.stress.aphi.vals(1)),']');
    vals(13)=hDV.data.stress.aphi.sigvals(1);
    hDV.data.distributions.deterministicVals{13} = hDV.data.stress.aphi.vals(1);
    deterministicValNums(13)= hDV.data.stress.aphi.vals(1);
    
    if hDV.data.stress.aphi.use == 11
        noshow = [2 3 4 9 11 12];
        % Reset variability for horizontal stresses and friction
        % coefficients to zero
        hDV.data.sigvals(2) = 0;
        hDV.data.sigvals(3) = 0;
    elseif hDV.data.stress.aphi.use == 12            
        noshow = [3 4 9 11 12] ;  % Modified A-phi model; show shmin
        % Reset variability for max horiz stress and ref mu to zero
        hDV.data.sigvals(3) = 0;    
    end
end

dataentryProbabilistic(hDV,title,txt,vals,noshow,deterministicValNums);
hDV.data.distributions.distirbutionTxt=txt; % save text


end


