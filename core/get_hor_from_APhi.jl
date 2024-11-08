module GetHorFromAPhi

export getHorFromAPhi

function getHorFromAPhi(APhi, mu, Sh_initial, Sv, pp0, min_hor_stress=nothing)
    println("Getting horizontal stresses using A-Phi model...")

    # Determine stress state (n) based on APhi
    n = if 0 <= APhi < 1
        0  # Normal faulting regime
    elseif 1 <= APhi < 2
        1  # Strike-slip regime
    elseif 2 <= APhi <= 3
        2  # Reverse faulting regime
    else
        error("APhi must be between 0-3. The APhi value you provided was $APhi")
    end

    # Calculate Phi based on APhi and n
    Phi = (APhi - (n + 0.5)) / (-1)^n + 0.5

    if min_hor_stress !== nothing
        println("Calculating horizontal stresses using provided minimum horizontal stress.")
        # Use min_hor_stress in calculations
        Sh = min_hor_stress
        # Effective stresses (assuming pore pressure is pp0)
        Sh_eff = Sh - pp0
        Sv_eff = Sv - pp0

        if n == 0
            # Normal faulting regime (Sv > SH > Sh)
            # SH_eff = Phi * (Sv_eff - Sh_eff) + Sh_eff
            SH_eff = Phi * (Sv_eff - Sh_eff) + Sh_eff
            SH = SH_eff + pp0
        elseif n == 1
            # Strike-slip regime (SH > Sv > Sh)
            # Sv_eff = Phi * (SH_eff - Sh_eff) + Sh_eff
            # Rearranged:
            # SH_eff = Sh_eff + (Sv_eff - Sh_eff) / Phi
            SH_eff = Sh_eff + (Sv_eff - Sh_eff) / Phi
            SH = SH_eff + pp0
        elseif n == 2
            # Reverse faulting regime (SH > Sh > Sv)
            # Sh_eff = Phi * (SH_eff - Sv_eff) + Sv_eff
            # Rearranged:
            # SH_eff = Sv_eff + (Sh_eff - Sv_eff) / Phi
            SH_eff = Sv_eff + (Sh_eff - Sv_eff) / Phi
            SH = SH_eff + pp0
        else
            error("Invalid faulting regime n = $n")
        end

    else
        println("Calculating horizontal stresses based on critical stress assumption.")
        # Use the original A-Phi model with critical stress assumption
        if mu > 0
            k = (mu + sqrt(1 + mu^2))^2
            Sv_eff = Sv - pp0  # Effective vertical stress

            if n == 0
                # Normal faulting regime
                Sh = (Sv - pp0) / k + pp0
                SH = Phi * (Sv - Sh) + Sh
            elseif n == 1
                # Strike-slip regime
                # Solve for Sh_eff and SH_eff
                denominator = Phi * (k - 1) + 1
                if denominator == 0
                    error("Denominator in Sh_eff calculation is zero")
                end
                Sh_eff = Sv_eff / denominator
                SH_eff = k * Sh_eff
                Sh = Sh_eff + pp0
                SH = SH_eff + pp0
            elseif n == 2
                # Reverse faulting regime
                Sh_eff = ([Phi * (k - 1) + 1]) * Sv_eff
                Sh = Sh_eff + pp0
                SH = k * Sv_eff + pp0
            else
                error("Invalid faulting regime n = $n")
            end
        else
            # In the case where mu <= 0, set SH and Sh equal to Sv
            SH = Sv
            Sh = Sv
        end
    end

    return SH, Sh  # Return the maximum and minimum horizontal stresses
end

end  # End of module
