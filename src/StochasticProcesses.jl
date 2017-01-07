module StochasticProcesses

include("results.jl")

export SimResult, CumsimResult
export initial, cumsim, distribution, sim
export BrownianMotion,
       BrownianMotionWithDrift,
       GeometricBrownianMotion,
       ItoIntegral,
       ItoProcess,
       SDE,
       CompositeProcess,
       AItoProcess

using Distributions

import Base.convert, Base.size, Base.rand

# Internal stochastic process interface
abstract AProcess

include("solvers.jl")

function Base.rand{P <: AProcess}(process::P, t, k::Int=1)
  # todo: better way to code this?
  if method_exists(distribution, (typeof(process), Float64))
    return Base.rand(distribution(process, t[end]), k)
  end
  return sim(process, t, k)
end

"Ito Process(f,y0): process satisfying equation dy = f(t, dt, b, db, y)."
immutable ItoProcess{F, Y} <: AProcess
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
immutable CompositeProcess{P <: AProcess, F} <: AProcess
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
abstract AItoProcess <: AProcess

init!{P <: AItoProcess}(process::P, state) = 
    init!(convert(ItoProcess, process), state)

next!{P <: AItoProcess}(process::P, state, i) = 
    next!(convert(ItoProcess, process), state, i)

size{P <: AItoProcess}(process::P) = size(convert(ItoProcess, process))

initial{P <: AItoProcess}(process::P) = initial(convert(ItoProcess, process))

# BrownianMotion

immutable BrownianMotion{Y} <: AItoProcess
  y0::Y
end

BrownianMotion() = BrownianMotion(0.0)

convert(::Type{ItoProcess}, bm::BrownianMotion) = 
    ItoProcess((t,dt,b,db,y)->db, bm.y0)

function distribution(bm::BrownianMotion, t) 
  if t == 0 
    Constant(bm.y0)
  elseif ndims(bm.y0) == 0
    Normal(bm.y0, sqrt(t))
  else 
    MvNormal(bm.y0, eye(length(bm.y0)) * t)
  end
end

# BrownianMotionWithDrift

immutable BrownianMotionWithDrift{Y,S} <: AItoProcess
  mu::Y
  sigma::S
  y0::Y
end

BrownianMotionWithDrift(mu, sigma) = BrownianMotionWithDrift(mu, sigma, 0.0)


convert(::Type{ItoProcess}, bm::BrownianMotionWithDrift) = 
    ItoProcess((t,dt,b,db,y)-> (bm.mu * dt) .+ (bm.sigma * db), bm.y0)

function distribution(bm::BrownianMotionWithDrift, t)
    if t == 0 
      Constant(bm.y0) 
    elseif ndims(bm.y0) == 0
      Normal(bm.y0 + bm.mu * t, bm.sigma * sqrt(t))
    else
      MvNormal(bm.y0 + bm.mu * t, bm.sigma *bm.sigma' * t)
    end
end

# Geometric Brownian Motion

type GeometricBrownianMotion <: AItoProcess
    mu::Float64
    sigma::Float64
    y0::Float64
end

convert(::Type{ItoProcess}, bm::GeometricBrownianMotion) = 
    ItoProcess((t,dt,b,db,y)-> (bm.mu * dt) * y + bm.sigma * (y .* db), bm.y0)

distribution(bm::GeometricBrownianMotion, t) = 
    t == 0 ? Constant(bm.y0) :
    LogNormal(log(bm.y0) + t * (bm.mu - bm.sigma^2 / 2), bm.sigma * sqrt(t))


# ItoIntegral

immutable ItoIntegral{F} <: AItoProcess
  f::F
end

convert{F}(::Type{ItoProcess}, ii::ItoIntegral{F}) = 
    ItoProcess((t,dt,b,db,y)->ii.f(t).*db, 0.)

immutable SDE{F,G} <: AItoProcess
  f::F
  g::G
  y0::Float64
end

convert{F,G}(::Type{ItoProcess}, sde::SDE{F,G}) = 
    ItoProcess((t,dt,b,db,y)->sde.f(t, y).*dt + sde.g(t, y).*db, sde.y0)


# Utilities

# Distribution with only one value
type Constant
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