% back door during setup to increase well or fault limit
% % include a .txt file in startup directory
% If the .txt file isn't in the same directory as the program, it won't be
% loaded
% running these 3 lines in matlab and putting this in the main directory
% edits the variables
% nwells_max=200; % 100 is default maximum
% NFAULTSMAX=600; % 500 faults is default maximum
% axsBgrnd = [1 , 1, 1] .*.4 % .8 is default RGB color
% mohrCircleDotSize = 4;  % 10 is default dot size
% save('./FSPStartup.mat','nwells_max','NFAULTSMAX','axsBgrnd','mohrCircleDotSize')


% 
function hDV = checkForStartupTxtFile(hDV)

c=dir('.');
b=struct2cell(c);% is this filename in current directory?
if any(ismember(b(1,:),'FSPStartupSpecifications.txt'))
% if mat file exists, load it
% load('./FSPStartup.mat')
fid3=fopen('./FSPStartupSpecifications.txt');

%     Read in as a string until you hit newline
out = textscan(fid3,'%s','delimiter','\n');
%     Get rid of empty lines
% out = cellfun(@(x) isempty(x),out{1});
out = out{1}(cellfun(@(x) ~isempty(x),out{1}));
fclose(fid3);
lines = length(out);
for i = 1:lines-1
    eval(out{i})
end
    
end

end








