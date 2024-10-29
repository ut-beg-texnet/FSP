module GetHorFromAPhi

export getHorFromAPhi

# Calculates horizontal stresses from critical stress assumptions and input of APhi and reference friction angle mu.
# Inputs: Aphi, mu, Sv (vertical stress), Sh (minmum horizontal stress), Pp (pore pressure), aphi_use (flag for model)
# Outputs: Horizontal total stresses SH and Sh


# Change this so it reads from the class object (data[:stress][:aphi][:use] is the value for aphi_use)
function getHorFromAPhi(APhi, mu, Sh, Sv, Pp, aphi_use)
    # determine stress state (n) based on APhi
    if 0 <= APhi < 1
        n = 0
    elseif 1 <= APhi < 2
        n = 1
    elseif 2 <= APhi <= 3
        n = 2
    else
        n = NaN
    end

    if !isnan(n)
        Phi = (APhi - (n + 0.5)) / (-1)^n + 0.5

        # if the 'aphi_use' flag is set to '12', use the modified A-Phi model
        if aphi_use == 12
            # Using Modified A-Phi model: Sh min as a model parameter instead of ref mu
            # effective stresses = total stresses - pore pressure
            Sh_eff = Sh - Pp # horizontal effective stress
            Sv_eff = Sv - Pp # vertical effective stress
            if n == 0
                SH = Phi * (Sv_eff - Sh_eff) + Sh_eff + Pp
            elseif n == 1
                SH = (Sv_eff - Sh_eff + Phi * Sh_eff) / Phi + Pp
            elseif n == 2
                SH = (Sh_eff - Sv_eff + Phi * Sv_eff) / Phi + Pp
            end
        else # aphi_use flag is not 12, use the original A-Phi model (based on mu)
            if mu > 0
                k = (mu + sqrt(1 + mu^2))^2
                if n == 0
                    Sh = (Sv - Pp) / k + Pp
                    SH = Phi * (Sv - Sh) + Sh
                elseif n == 1
                    A = [1.0 -k; Phi (1 - Phi)]
                    b = [Pp - k * Pp; Sv]
                    x = A \ b
                    SH = x[1]
                    Sh = x[2]
                elseif n == 2
                    SH = k * (Sv - Pp) + Pp
                    Sh = Phi * (SH - Sv) + Sv
                end
            else
                SH = Sv
                Sh = Sv
            end # mu > 0

        end # aphi_use

    else
        error("Check Aphi [0,3] range. You have APhi = $APhi")
    end

    return SH, Sh # returns SH (maximum horizontal stress) and Sh (minimum horizontal stress)
end

end # module
