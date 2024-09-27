% callback for when a well is selected in dropdown menu of model inputs
    function callwellhighlight2(src,~,hDV)
        valDropdown = get(src,'Value') ;
        numOtherOptions=2; % 1 ia all, 2 is none
        if  valDropdown==1 % then All was selected
            set(hDV.plotdata.pflot.pl(:),'linewidth',2)
            set(hDV.plotdata.pflot.htx(:),'fontsize',16,'fontweight','normal','visible','on')
        elseif  valDropdown==2 % then None was selected
            set(hDV.plotdata.pflot.pl(:),'linewidth',2)
            set(hDV.plotdata.pflot.htx(:),'fontsize',16,'fontweight','normal','visible','off')
            
        else % one well was selected
            set(hDV.plotdata.pflot.pl(:),'linewidth',2);
            % make bold line width corresponding to that well
            set(hDV.plotdata.pflot.pl(valDropdown-numOtherOptions),'linewidth',5);
            % bring line to top
            uistack(hDV.plotdata.pflot.pl(valDropdown-numOtherOptions),'top')
            % make bold font text
            set(hDV.plotdata.pflot.htx(:),'fontsize',16,'fontweight','normal')
            set(hDV.plotdata.pflot.htx(valDropdown-numOtherOptions),'fontsize',16+5,'fontweight','bold','visible','on')
            uistack(hDV.plotdata.pflot.htx(valDropdown-numOtherOptions),'top')% bring to top
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%