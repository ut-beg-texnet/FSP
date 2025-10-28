function [hDV,flagHydrology]=checkHydrologyData(hDV)

% check read in data to make sure it's compatible, and numbers/strings
% where it should be.
flagHydrology=0;
stringsToCheck=hDV.data.reservoir.stringsImportedHydrology;
[numrows,numcol]=size(stringsToCheck);
columnIsNumber=[1,1,1,1]; % x,y,p,t, should all be numbers

if numcol~=4
    errorWindow1=errordlg(cat(2,'you loaded in ',num2str(numcol),' columns of data, but there should be 4'));
    centerFigure(hDV.hfig,errorWindow1);
    flagHydrology=1;
end

errorcount=0;
% cycle over each row and column
for r = 1:numrows
    for c=1:numcol
        %        r = callbackdata.Indices(1); % row edited
        %         c = callbackdata.Indices(2); % column edited
        enteredVal=stringsToCheck{r,c};
        [~, numericFlag] = str2num(enteredVal); % try and convert to number, if unsuccessful, then numericFlag=0
        %         d=get(hObject,'Data');
        if numericFlag && columnIsNumber(c)  % you entered numbers and column takes numbers
            % do nothing, looks good
        elseif  ~ columnIsNumber(c) % you entered characters or numbers and column takes either
            % do nothing, looks good
            % I don't care if well name is numbers or letters or both
        elseif ~numericFlag % error
            errorWindow1= errordlg(cat(2,'in row ',num2str(r),' and column ',num2str(c),' you entered: ',enteredVal,' but expected a number'));
            centerFigure(hDV.hfig,errorWindow1);
            flagHydrology=1;
            errorcount=errorcount+1;
        else
            msgWindow1=msgbox(['something went wrong in checking entered value: ',enteredVal], 'entered value error','error'); % this should never be executed
            centerFigure(hDV.hfig,msgWindow1);
             flagHydrology=1;
        end
        
        if errorcount>10
             errorWindow2= errordlg(cat(2,'too many errors in hydology data, so couldn''t load it. Try cleaning it up. '));
            centerFigure(hDV.hfig,errorWindow2);
            return
        end
        
    end
end

%                 [rowsNaN,colsNaN]=find(isnan(hDV.data.fault.file));
%                 allOneString = sprintf('%.0f,' , rowsNaN);
%                 allOneString2 = sprintf('%.0f,' , colsNaN);
%                 errorWindow1= errordlg(cat(2,'you entered something that is Not a Number in row(s) ',allOneString,' and columns ',allOneString2,...
%                     ' (respectively). Fix this manually or re-load the csv file before proceeding.'));
%                 centerFigure(hfig,errorWindow1);


end