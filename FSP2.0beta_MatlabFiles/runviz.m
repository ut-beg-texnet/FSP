% D. Pais, R. Walsh and S.P. Lele
% ExxonMobil Upstream Research Company 2016+
% Stanford University 2016+

% Feb 11 2016 (v0.1)
% Apr 12 2016 (v1.0 legacy)

% August 10 2016 (v0.7)
% August 23 2016 (v0.8) 
% August 25 2016 (v0.81)
% Sept 15 2016 (v0.9) beta test compile for 5 SCITS affiliates
% Nov 2016 (v0.96) 
% Feb 2017 (V0.99)
% V 1.03 by Rall outside of Stanford thanks to apache 2.0 license
% V 1.05 change aphi treatment
% V 1.06 save .fig image 
% V 1.07 fix 0 inj rate and 0 year bugs, re-include black command window in
% deployed version by making it a console application in .prj settings

% Integrated Visualization for Disposal Well Induced Seismicity
% Quantitative Risk Analysis

model = 'FSP 2.0';
%add the path only on deployed application
if ~isdeployed
    addpath(genpath(pwd));
    vers = ['Coding Version ',model];
else
    vers = model ;
end


hSV = surfaceviz(vers) ; 



% toolbox dependency notes: 
% Statistics and machine learning toolbox: 
% ECDF function
% 
% Neural Network toolbox
% Combvec function
%
% image processing toolobox 
% imshow
%










