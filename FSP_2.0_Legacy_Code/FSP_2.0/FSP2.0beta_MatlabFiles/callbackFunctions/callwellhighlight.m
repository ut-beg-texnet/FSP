% callback for when a well is selected in dropdown menu of model inputs
    function callwellhighlight(src,~,hDV)
        valDropdown = get(src,'Value') ;
               numOtherOptions=2; % 1 ia all, 2 is none
        if  valDropdown==1 % then All was selected
            set(hDV.plotdata.inputMap.hwlocs2(:),'linewidth',2)
            set(hDV.plotdata.inputMap.hwellLabel(:),'fontsize',20,'fontweight','normal','visible','on')
            
       elseif  valDropdown==2 % then None was selected
            set(hDV.plotdata.inputMap.hwlocs2(:),'linewidth',2)
            set(hDV.plotdata.inputMap.hwellLabel(:),'fontsize',20,'fontweight','normal','visible','off')
            
        else % one well was selected
            % make bold line width corresponding to that well
            set(hDV.plotdata.inputMap.hwlocs2(valDropdown-numOtherOptions),'linewidth',5);
            notSelectedWell=[1:1:hDV.data.nwells]; % find wells not selected
            notSelectedWell(notSelectedWell==valDropdown-numOtherOptions)=[];
            set(hDV.plotdata.inputMap.hwlocs2(notSelectedWell),'linewidth',2) % make sure they have normal line width
            % bring line to top
            uistack(hDV.plotdata.inputMap.hwlocs2(valDropdown-numOtherOptions),'top')
            % make bold font text
            set(hDV.plotdata.inputMap.hwellLabel(valDropdown-numOtherOptions),'fontsize',20+5,'fontweight','bold','visible','on')
            set(hDV.plotdata.inputMap.hwellLabel(notSelectedWell),'fontsize',20,'fontweight','normal')
            uistack(hDV.plotdata.inputMap.hwellLabel(valDropdown-numOtherOptions),'top')% bring to top
            
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% callback %%%%%%%%%%%%%%%%