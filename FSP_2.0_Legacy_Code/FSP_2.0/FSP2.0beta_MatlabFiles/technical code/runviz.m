% D. Pais & R. Walsh
% ExxonMobil Upstream Research Company 2016+
% Stanford University 2016+

% Feb 11 2016 (v0.1)
% Apr 12 2016 (v1.0 legacy)

% August 10 2016 (v0.7)
% August 23 2016 (v0.8) 
% August 25 2016 (v0.81)
% Sept 15 2016 (v0.9) beta test compile for 5 SCITS affiliates
% Nov 2016 (v0.96) 

% Integrated Visualization for Disposal Well Induced Seismicity
% Quantitative Risk Analysis

%add the path only on deployed application
if ~isdeployed
    addpath(genpath(pwd));
    vers = 'Running Code';
else
    vers = '0.98.8 (Beta)' ;
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












