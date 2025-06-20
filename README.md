## BackwardDifferenceFormula.jl

A flexible and efficient tool to evaluate Backward Difference Formula (BDF) coefficients for variable time ste and arbitrary order.

The function `bfd_coefficients()` can be used to determine the BDF coefficient provided the input `tshift` arrays is provided. This array contains the value of relative time from the oldest point to the newest. Since it's a relative values, the newest time has value `0`. `tshift` can be either a `SVector` or a `Tuple`.

### Quick start

In package mode (type `]` in Julia's REPL):
1) `add https://github.com/tduretz/BackwardDifferenceFormula.jl`

Go back to evaluation mode (press backspace)

2) `using BackwardDifferenceFormula` 

For the best experience, use [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl), as in the examples below, or `Tuple`s for the `tshift` argument.

#### BDF 1

```julia 
tshift = @SVector [-1.0, 0.0] # assumes Δt = 1
    
coeff = bdf_coefficients(tshift)

display(coeff)
```
The resulting coefficients should be those of backward-Euler, `1/Δt` and `-1/Δt`.

#### BDF 2

```julia 
tshift = @SVector [-2.0, -1.0, 0.0] # assumes Δt = 1 and constant
    
coeff = bdf_coefficients(tshift)

display(coeff)
```
The resulting coefficients should be, for `Δt = 1` , `3/2`, `-2`, `1/2`.

### Example: Visco-Elasto-Plastic rheology integration

This example depicts the integration of visco-elasto-plastic rheology (Maxwell assembly) using various order BDF. This examples uses constant step. The code is provided [here](examples/viscoelastoplastic.jl).

The solution is a stress-time curve which looks like this:
![](/images/viscoelastoplastic_run.svg).

The theoretical rates of convergence of various order BDF are achieved:
![](/images/viscoelastoplastic.svg).


