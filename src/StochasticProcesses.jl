module StochasticProcesses

export cumsim, sim, initial
export SimResult, CumsimResult

export AItoProcess,
       BrownianBridge,
       BrownianMotion,
       BrownianMotionWithDrift,
       CompositeProcess,
       FirstTime,
       GeometricBrownianMotion,
       ItoIntegral,
       ItoProcess,
       PoissonProcess,
       SDE

# Basic processes information
export distribution, solution

# Internal stochastic process interface
abstract type AProcess end

include("brownian_bridge.jl")
include("first_time.jl")
include("poisson.jl")
include("results.jl")
include("solvers.jl")

using Distributions
import Base.convert, Base.size, Base.rand

function Base.rand{P <: AProcess}(process::P, t, k::Int=1)
  # todo: better way to code this?
  if method_exists(distribution, (typeof(process), Float64))
    return Base.rand(distribution(process, t[end]), k)
  end
  return sim(process, t, k)
end

"Ito Process(f,y0): process satisfying equation dy = f(t, dt, b, db, y)."
struct ItoProcess{F, Y} <: AProcess
  f::F
  y0::Y
end

function init!{F,Y}(process::ItoProcess{F,Y}, g::Generator)
  g.y = repeat(vwrap(initial(process)), outer=(1, g.k))
  g.y
end

@inline function next!{F,Y}(process::ItoProcess{F,Y}, g::Generator, i)
  dy = process.f(g.t[i-1], g.dt[i-1], g.b, g.db, g.y)
  add!(g.y, dy)
  g.y
end

initial(process::ItoProcess) = process.y0

# CompositeProcess
# The resulting process is f(t, y), where y is the value of p
struct CompositeProcess{P <: AProcess, F} <: AProcess
  p::P
  f::F
end

init!{P, F}(process::CompositeProcess{P,F}, g::Generator) = 
    process.f.(0., init!(process.p, g))

@inline function next!{P, F}(process::CompositeProcess{P,F}, g::Generator, i) 
  next!(process.p, g, i)
  process.f.(g.t[i], g.y)
end

initial{P, F}(process::CompositeProcess{P,F}) = initial(process.p)

size{P, F}(process::CompositeProcess{P,F}) = size(process.p)

"A process that is convertible to ItoProcess"
abstract type AItoProcess <: AProcess end

init!{P <: AItoProcess}(process::P, state) = 
    init!(convert(ItoProcess, process), state)

next!{P <: AItoProcess}(process::P, state, i) = 
    next!(convert(ItoProcess, process), state, i)

size{P <: AItoProcess}(process::P) = size(convert(ItoProcess, process))

initial{P <: AItoProcess}(process::P) = initial(convert(ItoProcess, process))

include("basic_processes.jl")

# Utilities

# Distribution with only one value
mutable struct Constant
  x
end

Base.rand(c::Constant, k) = fill(c.x, k)

function collapse_but_first(A::AbstractArray)
    s = size(A)
    dims = tuple(find([d == 1 && i != 1 for (i,d) in enumerate(s)])...)
    return squeeze(A,dims)
end

vwrap(x::AbstractArray) = x
vwrap(x::Number) = [x]

@inline function add!(y, dy)
  for i in eachindex(y)
    @inbounds y[i] += dy[i]
  end
end

@inline function cmul!(y, dy)
  for i in eachindex(y)
    @inbounds y[i] *= dy
  end
end


end