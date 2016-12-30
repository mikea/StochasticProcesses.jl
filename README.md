# StochasticProcesses.jl

A Julia package for (continuous) stochastic processes.

## Examples

```julia

julia> using StochasticProcesses;

julia> cumsim(BrownianMotion(), linspace(0, 1, 5))
5×1 Array{Float64,2}:
  0.0      
 -0.705048 
 -0.0780865
 -0.615589 
 -0.93188  

julia> cumsim(BrownianMotion(), linspace(0, 1, 5), 3)
5×3 Array{Float64,2}:
  0.0        0.0        0.0      
 -0.815025  -0.331914  -0.146447 
 -0.643032   0.230647  -0.528345 
 -1.04061   -0.713817  -0.0711855
 -0.226368  -0.417196  -0.203268 

julia> sim(BrownianMotionWithDrift(10, 10, 100), linspace(0, 1, 1000), 3)
3-element Array{Float64,1}:
  97.8373
 107.271 
 105.589 

julia> distribution(GeometricBrownianMotion(.1, .3, 100), 10)
Distributions.LogNormal{Float64}(μ=5.155170185988092, σ=0.9486832980505138)
```

`rand(process, t, k==1)` - returns `k` samples of `process` at the end of the time grid `t`.
Uses precise analytical distribution if available, or `sim`.

```julia
julia> rand(BrownianMotionWithDrift(10, 10, 100), linspace(0, 1, 1000), 3)
3-element Array{Float64,1}:
 115.006
 108.605
 101.62 
```
