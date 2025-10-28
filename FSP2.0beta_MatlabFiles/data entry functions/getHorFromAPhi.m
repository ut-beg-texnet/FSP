function [SH,Sh] = getHorFromAPhi(APhi,mu,Sh,Sv,Pp,aphi_use)

% Darren Pais
% ExxonMobil URC 2016
% Function calculates horizontal stresses from critical stress assumptions
% and input of APhi and reference friction angle mu. Total vertical stress
% Sv and pore pressure Pp are the other inputs. 
%
% "Modified A-Phi" calculations added by Suvrat Lele  March 2018

% Outputs are the horizontal total stresses SH and Sh
% Look at Hurd and Zoback 2012 as a reference


%calculate n and Phi from APhi
if APhi>=0 && APhi<1
    n=0;
elseif APhi>=1 && APhi<2
    n=1;
elseif APhi>=2 && APhi<=3
    n=2;
else
    n=NaN;
end


if ~isnan(n)
    
    Phi = (APhi-(n+0.5))/(-1)^n+0.5 ;
    
    if aphi_use==12
        
        % Using Modified A-Phi model: 
        % use shmin as a model parameter instead of ref mu
        %
        Sh_eff = Sh-Pp;
        Sv_eff = Sv-Pp;
        %
        switch n
            case 0
                SH = Phi*(Sv_eff-Sh_eff)+ Sh_eff + Pp;
            case 1
                SH = (Sv_eff-Sh_eff+Phi*Sh_eff)/Phi + Pp;
            case 2
                SH = (Sh_eff-Sv_eff+Phi*Sv_eff)/Phi + Pp;
        end
        
    else
        
        if mu > 0
            k   = (mu + sqrt(1+mu.^2))^2;
            switch n
                case 0
                    Sh = (Sv-Pp)/k + Pp;
                    SH = Phi*(Sv-Sh)+Sh;
                case 1
                    A=[1 -k ; Phi (1-Phi)]; b=[Pp-k*Pp ; Sv]; x=A\b ;
                    SH=x(1) ; Sh=x(2) ;
                case 2
                    SH = k*(Sv-Pp)+Pp ;
                    Sh = Phi*(SH-Sv)+Sv ;
            end
        else
            SH = Sv;
            Sh = Sv;
        end % mu >0
        
    end % aphi_use
    
else
    edlgBox= errordlg(['Check Aphi [0,3] range. You have APhi =',num2str(APhi)]) ;
    try
        centerFigure(hDV.hfig,edlgBox);
    catch
        centerFigure(gcf,edlgBox);
    end
end


end