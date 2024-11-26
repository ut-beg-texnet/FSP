module GetHorFromAPhi

export getHorFromAPhi

function getHorFromAPhi(APhi, mu, Sv, pp0, min_hor_stress=nothing)
    println("getHorFromAPhi inputs: APhi=$APhi, mu=$mu, Sv=$Sv, pp0=$pp0, min_hor_stress=$min_hor_stress")

    # Determine stress state (n) based on APhi
    if 0 <= APhi && APhi < 1
        n = 0  # Normal faulting regime
    elseif 1 <= APhi && APhi < 2
        n = 1  # Strike-slip regime
    elseif 2 <= APhi && APhi <= 3
        n = 2  # Reverse faulting regime
    else
        error("APhi must be between 0-3. The APhi value you provided was $APhi")
    end

    # Calculate Phi based on APhi and n
    Phi = (APhi - (n + 0.5)) / (-1)^n + 0.5 # VERIFY PARENTHESES
    #Phi = (0.58 - (0 + 0.5)) / (-1)^0 + 0.5


    if min_hor_stress !== nothing
        println("Calculating horizontal stresses using provided minimum horizontal stress.")
        # Use min_hor_stress in calculations
        Sh = min_hor_stress
        println("Sh = $Sh")
        # Effective stresses
        Sh_eff = Sh - pp0
        Sv_eff = Sv - pp0


        if n == 0
            # Normal faulting regime (Sv > SH > Sh)
            #SH_eff = Phi * (Sv_eff - Sh_eff) + Sh_eff
            #SH = SH_eff + pp0
            SH = Phi*(Sv_eff-Sh_eff)+ Sh_eff + pp0
        elseif n == 1
            # Strike-slip regime (SH > Sv > Sh)
            #SH_eff = Sh_eff + (Sv_eff - Sh_eff) / Phi
            #SH_eff = Sh_eff + (Sv_eff - Sh_eff) / Phi
            #SH = SH_eff + pp0
            SH = (Sv_eff-Sh_eff+Phi*Sh_eff)/Phi + pp0
        elseif n == 2
            # Reverse faulting regime (SH > Sh > Sv)
            #SH_eff = Sv_eff + (Sh_eff - Sv_eff) / Phi
            #SH = SH_eff + pp0
            SH = (Sh_eff-Sv_eff+Phi*Sv_eff)/Phi + pp0
        else
            error("Invalid faulting regime n = $n")
        end

    else
        println("Calculating horizontal stresses based on critical stress assumption.")
        if mu > 0
            k = (mu + sqrt(1 + mu^2))^2
            #Sv_eff = Sv - pp0  # Effective vertical stress

            if n == 0
                # Normal faulting regime
                #Sh_eff = (Sv_eff) / k
                #Sh = Sh_eff + pp0
                #SH = Phi * (Sv - Sh) + Sh
                Sh = (Sv-pp0)/k + pp0
                SH = Phi*(Sv-Sh)+Sh
            elseif n == 1
                # Strike-slip regime
                # Solve the linear system A * x = b
                A = [1.0 -k; Phi (1 - Phi)]
                b = [pp0 - k*pp0 ; Sv]
                x = A \ b
                SH = x[1]
                Sh = x[2]
            elseif n == 2
                # Reverse faulting regime
                #H_eff = k * Sv_eff
                #SH = SH_eff + pp0
                #Sh = Phi * (SH - Sv) + Sv
                SH = k*(Sv-pp0) + pp0
                Sh = Phi*(SH-Sv)+Sv
            else
                error("Invalid faulting regime n = $n")
            end
        else
            # In the case where mu <= 0, set SH and Sh equal to Sv
            SH = Sv
            Sh = Sv
            
        end
    end

    println("End of getHorFromAPhi function.")

    return SH, Sh  # Return the maximum and minimum horizontal stresses
end

end  # End of module
