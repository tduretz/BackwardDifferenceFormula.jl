using Test, StaticArrays
using BackwardDifferenceFormula

function main()

    @testset "BDF1    " begin

        tshift = @SVector([-1.0, 0.0])
    
        coeff = bdf_coefficients(tshift)
        @show coeff
        # box = BBox(o2, 2, 4)
        # @test inside(p, box)
    end
    
    return
end

main()