# StochasticProcesses.jl

A Julia package for (continuous) stochastic processes. Monte-Carlo simulation and basic facts.

Demo notebooks: 
* [examples/sde.ipynb](examples/sde.ipynb) - demostration of basic stochastic differential equations properties.
* [examples/useful_processes.ipynb)(examples/useful_process.ipynb) - demonstration of basic processes and their properties.
* [examples/black-scholes.ipynb](examples/black-scholes.ipynb) - Black Scholes formula demonstration.

## Overview

```julia

julia> using StochasticProcesses;

```

### Basic Processes

`BrownianMotion(y0=0)`, `BrownianMotionWithDrift(mu, sigma, y0)`, `GeometricBrownianMotion(mu, sigma, y0)`.


### Ito Processes

`ItoIntegral(f)` - solution of `dy = f(y) * dB` with `y(0) = 0`.

`SDE(f, g, y0)` - solution of `dy = f(t, y) * dt + g(t, y) * dB`.

`CompositeProcess(process, f)` - `f(t, y)` applied to the `process`.

### Simulating Paths

`cumsim(process, t, k==1)` returns `k` paths of `process` simulation on the time grid `t`.

```julia
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
```

### Simulated Sampling

`sim(process, t, k)` returns `k` end results of `process` simulation on the time grid `t.`
Faster and more memory efficient that `cumsim.`

```julia
julia> sim(BrownianMotionWithDrift(10, 10, 100), linspace(0, 1, 1000), 3)
3-element Array{Float64,1}:
  97.8373
 107.271 
 105.589 
 ```
 
### Distribution
 
`distribution(process, t)` returns probability distribution of `process` value at time `t.`
Only available for: `BrownianMotion,` `BrownianMotionWithDrift,` `GeometricBrownianMotion.`

 ```julia
julia> distribution(GeometricBrownianMotion(.1, .3, 100), 10)
Distributions.LogNormal{Float64}(μ=5.155170185988092, σ=0.9486832980505138)
```

### Efficient Sampling

`rand(process, t, k==1)` returns `k` samples of `process` at the end of the time grid `t.`
Uses precise analytical distribution if available, or `sim.`

```julia
julia> rand(BrownianMotionWithDrift(10, 10, 100), linspace(0, 1, 1000), 3)
3-element Array{Float64,1}:
 115.006
 108.605
 101.62 
```
