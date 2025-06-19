using Plots, StaticArrays, BackwardDifferenceFormula, LinearAlgebra

function viscoelastoplastic(bdf, n)
    
    η      = 10.0
    G      = 1.0
    Δt_ref = 1.0/n
    Nt     = 100*n
    ε̇      = 1.0
    C      = 10

    t0    = 0.
    τ     = 2η*ε̇*(1 .- exp.(-G/η.*t0))
    τ0    = 2η*ε̇*(1 .- exp.(-G/η.*(t0-1*Δt_ref)))
    τ00   = 2η*ε̇*(1 .- exp.(-G/η.*(t0-2*Δt_ref)))
    τ000  = 2η*ε̇*(1 .- exp.(-G/η.*(t0-3*Δt_ref)))
    τ0000 = 2η*ε̇*(1 .- exp.(-G/η.*(t0-4*Δt_ref)))

    τvec = zeros(Nt)

    Δt = Δt_ref
    Δt0, Δt00, Δt000 = Δt, Δt, Δt_ref

    # Time stepping
    for it=1:Nt

        τ0000 = τ000
        τ000  = τ00
        τ00   = τ0
        τ0    = τ

        Δt000 = Δt00
        Δt00  = Δt0
        Δt0   = Δt

        Δt    = Δt_ref
        
        if bdf==1 
            tshift = @SVector([-Δt, 0.0])
            coeff  = bdf_coefficients(tshift)
            a, b, c, d, e = coeff[1], coeff[2], 0, 0, 0
        end

        if bdf==2 
            tshift = @SVector([-Δt-Δt0, -Δt, 0.0])
            coeff  = bdf_coefficients(tshift)
            a, b, c, d, e = coeff[1], coeff[2], coeff[3], 0, 0
        end

        if bdf==3 
            tshift = @SVector([-Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0.0])
            coeff  = bdf_coefficients(tshift)
            a, b, c, d, e = coeff[1], coeff[2], coeff[3], coeff[4], 0
        end

        if bdf==4 
            tshift = @SVector([-Δt-Δt0-Δt00-Δt000, -Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0.0])
            coeff  = bdf_coefficients(tshift)
            a, b, c, d, e = coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]
        end

        ηve   = ( 1/η + a/G )^-1
        τ     = 2*ηve* (ε̇ - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        F     = abs(τ) - C
        if F>0
            λ̇     = 0.5*F/ηve
            ε̇pl   = λ̇
            τ     = 2*ηve* (ε̇ - ε̇pl - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        end
        τvec[it] = τ
    end

    t   = (1:Nt).*Δt
    τvec_ana = 2η*ε̇*(1 .- exp.(-G/η.*t))
    τvec_ana[τvec_ana.>C] .= 10
    
    p = plot()
    p = plot!(t, τvec)
    p = plot!(t, τvec_ana )
    display(p)

    err = norm(τvec .- τvec_ana)/norm(τvec_ana)
    return err

end

let

res = zeros(4,4)

@info "BDF1"
res[1,1] = viscoelastoplastic(1,1)
res[1,2] = viscoelastoplastic(1,2)
res[1,3] = viscoelastoplastic(1,4)
res[1,4] = viscoelastoplastic(1,8)
@show res[1,:]

@info "BDF2"
res[2,1] = viscoelastoplastic(2,1)
res[2,2] = viscoelastoplastic(2,2)
res[2,3] = viscoelastoplastic(2,4)
res[2,4] = viscoelastoplastic(2,8)
@show res[2,:]

@info "BDF3"
res[3,1] = viscoelastoplastic(3,1)
res[3,2] = viscoelastoplastic(3,2)
res[3,3] = viscoelastoplastic(3,4)
res[3,4] = viscoelastoplastic(3,8)
@show res[3,:]

@info "BDF4"
res[4,1] = viscoelastoplastic(4,1)
res[4,2] = viscoelastoplastic(4,2)
res[4,3] = viscoelastoplastic(4,4)
res[4,4] = viscoelastoplastic(4,8)
@show res[4,:]

x = 1 ./[1 2 4 8]
p=plot(xlabel = "log(1/Δt)", ylabel="Stress error")
p=plot!(log10.(1 ./x[:]), log10.(res[1,:]), marker=:dot, label="BDF1")
p=plot!(log10.(1 ./x[:]), log10.(res[2,:]), marker=:dot, label="BDF2")
p=plot!(log10.(1 ./x[:]), log10.(res[3,:]), marker=:dot, label="BDF3")
p=plot!(log10.(1 ./x[:]), log10.(res[4,:]), marker=:dot, label="BDF4")
end