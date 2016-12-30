module StochasticProcesses

export cumsim, distribution, sim
export BrownianMotion,
       BrownianMotionWithDrift,
       GeometricBrownianMotion,
       ItoIntegral,
       ItoProcess,
       SDE,
       CompositeProcess,
       AItoProcess

using Distributions

import Base.convert

# Internal stochastic process implementation

abstract AProcess

type Generator
  # specified
  t
  n::Int64
  k::Int64

  # computed
  y::Vector{Float64}
  dt::Vector{Float64}
  sdt::Vector{Float64}
  b::Vector{Float64}
  db::Vector{Float64}

  function Generator(t, k::Int)
    g      = new(t, length(t), k)
    g.dt   = diff(t)
    g.sdt  = sqrt(g.dt)
    g.b    = zeros(k)
    g.db   = zeros(k)
    g.y    = zeros(k)
    return g
  end
end

@inline function add!(y, dy, k)
  for i in 1:k
    @inbounds y[i] += dy[i]
  end
end

@inline function next{P <: AProcess}(process::P, g::Generator, i::Int64)
  if i == 1
    return init!(process, g)
  end

  # Important to use previous db for process computation not to have bias.
  add!(g.b, g.db, g.k)
  randn!(g.db)
  for j in 1:g.k
     @inbounds g.db[j] *= g.sdt[i-1]
  end

  return next!(process, g, i)
end

function cumsim{P <: AProcess}(process::P, t, k::Int=1)
  gen = Generator(vec(t), k)
  result = zeros(gen.n, k)
  for i in 1:gen.n
    y = next(process, gen, i)
    for j in 1:k
      result[i,j] = y[j]
    end 
  end
  result
end

function sim{P <: AProcess}(process::P, t, k::Int=1)
  gen = Generator(t, k)
  for i in 1:(gen.n-1)
    next(process, gen, i)
  end
  next(process, gen, gen.n)
end

function Base.rand{P <: AProcess}(process::P, t, k::Int=1)
  # todo: better way to code this?
  if method_exists(distribution, (typeof(process), Float64))
    return Base.rand(distribution(process, t[end]), k)
  end
  return sim(process, t, k)
end

"Ito Process: dy = f(t, dt, b, db, y)."
immutable ItoProcess{F} <: AProcess
  f::F
  y0::Float64
end

function init!{F}(process::ItoProcess{F}, g::Generator)
  fill!(g.y, process.y0)
  g.y
end

@inline function next!{F}(process::ItoProcess{F}, g::Generator, i)
  dy = process.f(g.t[i-1], g.dt[i-1], g.b, g.db, g.y)
  add!(g.y, dy, g.k)
  g.y
end


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

# A process that is convertible to ItoProcess

abstract AItoProcess <: AProcess

init!{P <: AItoProcess}(process::P, state) = 
    init!(convert(ItoProcess, process), state)

next!{P <: AItoProcess}(process::P, state, i) = 
    next!(convert(ItoProcess, process), state, i)

# BrownianMotion

immutable BrownianMotion <: AItoProcess
  y0::Float64

  BrownianMotion() = new(0.0)
  BrownianMotion(y0) = new(y0)
end

convert(::Type{ItoProcess}, bm::BrownianMotion) = 
    ItoProcess((t,dt,b,db,y)->db, bm.y0)

distribution(bm::BrownianMotion, t) = t == 0 ? Constant(bm.y0) : Normal(bm.y0, sqrt(t))

# BrownianMotionWithDrift

immutable BrownianMotionWithDrift <: AItoProcess
  mu::Float64
  sigma::Float64
  y0::Float64

  BrownianMotionWithDrift(mu, sigma) = new(mu, sigma, 0.0)
  BrownianMotionWithDrift(mu, sigma, y0) = new(mu, sigma, y0)
end

convert(::Type{ItoProcess}, bm::BrownianMotionWithDrift) = 
    ItoProcess((t,dt,b,db,y)-> bm.mu * dt + bm.sigma * db, bm.y0)

distribution(bm::BrownianMotionWithDrift, t) = 
    t == 0 ? Constant(bm.y0) : Normal(bm.y0 + bm.mu * t, bm.sigma * sqrt(t))


# Geometric Brownian Motion

type GeometricBrownianMotion <: AItoProcess
    mu::Float64
    sigma::Float64
    y0::Float64
end

convert(::Type{ItoProcess}, bm::GeometricBrownianMotion) = 
    ItoProcess((t,dt,b,db,y)-> (bm.mu * dt) * y + bm.sigma * (y .* db), bm.y0)

distribution(bm::GeometricBrownianMotion, t) = 
    t == 0 ? Constant(0) :
    LogNormal(log(bm.y0) + t * (bm.mu - bm.sigma^2 / 2), bm.sigma * sqrt(t))


# ItoIntegral

immutable ItoIntegral <: AItoProcess
  f::Function
end

convert(::Type{ItoProcess}, ii::ItoIntegral) = 
    ItoProcess((t,dt,b,db,y)->ii.f(t).*db, 0.)

immutable SDE <: AItoProcess
  f::Function
  g::Function
  y0::Float64
end

convert(::Type{ItoProcess}, sde::SDE) = 
    ItoProcess((t,dt,b,db,y)->sde.f(t, y).*dt + sde.g(t, y).*db, sde.y0)


# Utilities

# Distribution with only one value
type Constant
  x
end

Base.rand(c::Constant, k) = Base.rand(c.x, k)

end