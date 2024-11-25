% change maximum X axis limit on hax, getting value from src


    function callmaxdpprob(src,~,hax,hDV)
        val = str2double(get(src,'string'));
        if val<=0 || isnan(val)
            try 
               val=max([floor( min(hDV.data.stress.vals(1:3)).*hDV.data.stress.vals(5)),1000]); 
            catch
            val=8000;
            end
        end
        xlim(hax,[0 val]) ; set(src,'string',num2str(val));        
    end
    
    
    