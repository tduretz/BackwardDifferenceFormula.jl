using Plots, StaticArrays, BackwardDifferenceFormula, LinearAlgebra

function viscoelastoplastic(bdf, n)

    # Physical parameters
    η      = 50.0
    G      = 1.0
    Δt_ref = 1.0/n
    Nt     = 100*n
    ε̇      = 1.0
    C      = 50

    # Initialise solution
    t0    = 0.
    τ     = 2η*ε̇*(1 .- exp.(-G/η.*t0))
    τ0    = 2η*ε̇*(1 .- exp.(-G/η.*(t0-1*Δt_ref)))
    τ00   = 2η*ε̇*(1 .- exp.(-G/η.*(t0-2*Δt_ref)))
    τ000  = 2η*ε̇*(1 .- exp.(-G/η.*(t0-3*Δt_ref)))
    τ0000 = 2η*ε̇*(1 .- exp.(-G/η.*(t0-4*Δt_ref)))
    τvec  = zeros(Nt)
    Δt, Δt0, Δt00, Δt000 = Δt_ref, Δt_ref, Δt_ref, Δt_ref

    # Select integrator
    if bdf == 1
        tshift = @SVector [-Δt, 0.0]
        coeff  = bdf_coefficients(tshift)
        a, b, c, d, e = coeff..., 0, 0, 0
    elseif bdf == 2
        tshift = @SVector [-Δt-Δt0, -Δt, 0.0]
        coeff  = bdf_coefficients(tshift)
        a, b, c, d, e = coeff..., 0, 0
    elseif bdf == 3
        tshift = @SVector [-Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0.0]
        coeff  = bdf_coefficients(tshift)
        a, b, c, d, e = coeff..., 0
    elseif bdf == 4
        tshift = @SVector [-Δt-Δt0-Δt00-Δt000, -Δt-Δt0-Δt00, -Δt-Δt0, -Δt, 0.0]
        coeff  = bdf_coefficients(tshift)
        a, b, c, d, e = coeff
    end

    # Time stepping
    for it in 1:Nt

        # Old stress values
        τ0000 = τ000
        τ000  = τ00
        τ00   = τ0
        τ0    = τ

        # Old time step values
        Δt000 = Δt00
        Δt00  = Δt0
        Δt0   = Δt

        # Visco-elastic trial
        ηve   = ( 1/η + a/G )^-1
        τ     = 2*ηve* (ε̇ - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        F     = abs(τ) - C

        # Visco-elastic-plastic corrector
        if F > 0
            λ̇     = 0.5*F/ηve
            ε̇pl   = λ̇
            τ     = 2*ηve* (ε̇ - ε̇pl - (b*τ0 + c*τ00 + d*τ000 + e*τ0000)/2/G )
        end
        τvec[it] = τ
    end

    # Exact solution
    t   = (1:Nt).*Δt
    τvec_ana = 2η*ε̇*(1 .- exp.(-G/η.*t))
    @views τvec_ana[τvec_ana.>C] .= C

    # Visualize
    p = plot(xlabel="Time", ylabel="Stress" )
    p = plot!(t, τvec, label="Numerics")
    p = plot!(t, τvec_ana, label="Exact" )
    display(p)

    return norm(τvec .- τvec_ana)/norm(τvec_ana)
end

let
    # Convergence test
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
    p=plot(xlabel = "log(1/Δt)", ylabel="Stress error", title="Convergence: visco-elasto-plasticity")
    p=plot!(log10.(1 ./x[:]), log10.(res[1,:]), marker=:dot, label="BDF1")
    p=plot!(log10.(1 ./x[:]), log10.(res[2,:]), marker=:dot, label="BDF2")
    p=plot!(log10.(1 ./x[:]), log10.(res[3,:]), marker=:dot, label="BDF3")
    p=plot!(log10.(1 ./x[:]), log10.(res[4,:]), marker=:dot, label="BDF4")
end
