using Test, StaticArrays
using BackwardDifferenceFormula

function main()

    @testset "BDF1    " begin

        tshift = @SVector([-1.0, 0.0])

        coeff = bdf_coefficients(tshift)

        # BDF1 (Backward Euler with Δt = 1.0)
        coeff_BDF1 = @SVector([1.0, -1.0])

        for i in eachindex(coeff)
            @test sum(coeff[i] .- coeff_BDF1[i]) ≈ 0
        end
    end

    @testset "BDF2    " begin

        tshift = @SVector([-2.0, -1.0, 0.0])

        coeff = bdf_coefficients(tshift)

        # BDF2 (with Δt = 1.0)
        coeff_BDF2 = @SVector([3/2, -4/3*3/2, 1/3*3/2])

        for i in eachindex(coeff)
            @test sum(coeff[i] .- coeff_BDF2[i]) ≈ 0
        end
    end

    @testset "BDF3    " begin

        tshift = @SVector([-3.0, -2.0, -1.0, 0.0])

        coeff = bdf_coefficients(tshift)

        # BDF3 (with Δt = 1.0)
        coeff_BDF3 = @SVector([11/6, -3, 3/2, -1/3])

        for i in eachindex(coeff)
            @test sum(coeff[i] .- coeff_BDF3[i]) ≈ 0 atol = 1e-15
        end
    end

    @testset "BDF4    " begin

        tshift = @SVector([-4.0, -3.0, -2.0, -1.0, 0.0])

        coeff = bdf_coefficients(tshift)

        # BDF4 (with Δt = 1.0)
        a = 25/12
        coeff_BDF4 = @SVector([25/12, -48/25*a, 36/25*a, -16/25*a, 3/25*a])

        for i in eachindex(coeff)
            @test sum(coeff[i] .- coeff_BDF4[i]) ≈ 0 atol = 1e-15
        end
    end

    return
end

main()
