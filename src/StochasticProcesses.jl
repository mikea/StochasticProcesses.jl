module StochasticProcesses

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

# Internal stochastic process implementation

abstract AProcess

type Generator
  # specified
  t::AbstractArray{Float64, 1}
  n::Int64
  "number of dimensions."
  k::Int64

  # computed
  d::Tuple
  y::Array{Float64}
  dt::Vector{Float64}
  sdt::Vector{Float64}
  b::Array{Float64}
  db::Array{Float64}

  function Generator(t, s, k::Int)
    g      = new(t, length(t), k)
    g.d    = (s..., k)
    g.dt   = diff(t)
    g.sdt  = sqrt(g.dt)
    g.b    = zeros(g.d)
    g.db   = zeros(g.d)
    g.y    = zeros(g.d)
    return g
  end
end

@inline function next{P <: AProcess}(process::P, g::Generator, i::Int64)
  if i == 1
    return init!(process, g)
  end

  # Important to use previous db for process computation not to have bias.
  add!(g.b, g.db)
  randn!(g.db)
  cmul!(g.db, g.sdt[i-1])
  next!(process, g, i)
end

function cumsim{P <: AProcess}(process::P, t, k::Int=1)
  s=size(vwrap(initial(process)))
  gen = Generator(vec(t), s, k)
  result = zeros(gen.n, s..., k)
  for i in 1:gen.n
    y = next(process, gen, i)
    result[i,:,:] .= y[:,:]
  end
  collapse_but_first(result)
end

function sim{P <: AProcess}(process::P, t, k::Int=1)
  s=size(vwrap(initial(process)))
  gen = Generator(t, s, k)
  for i in 1:(gen.n-1)
    next(process, gen, i)
  end
  collapse_but_first(next(process, gen, gen.n))
end

function Base.rand{P <: AProcess}(process::P, t, k::Int=1)
  # todo: better way to code this?
  if method_exists(distribution, (typeof(process), Float64))
    return Base.rand(distribution(process, t[end]), k)
  end
  return sim(process, t, k)
end

"Ito Process(f,y0): process satisfying equation dy = f(t, dt, b, db, y)."
immutable ItoProcess{F} <: AProcess
  f::F
  y0
end

function init!{F}(process::ItoProcess{F}, g::Generator)
  g.y = repeat(vwrap(initial(process)), outer=(1, g.k))
  g.y
end

@inline function next!{F}(process::ItoProcess{F}, g::Generator, i)
  dy = process.f(g.t[i-1], g.dt[i-1], g.b, g.db, g.y)
  add!(g.y, dy)
  g.y
end

initial{F}(process::ItoProcess{F}) = process.y0

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

immutable BrownianMotion <: AItoProcess
  y0

  BrownianMotion() = new(0.0)
  BrownianMotion(y0) = new(y0)
end

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

immutable BrownianMotionWithDrift <: AItoProcess
  mu
  sigma
  y0

  BrownianMotionWithDrift(mu, sigma) = new(mu, sigma, 0.0)
  BrownianMotionWithDrift(mu, sigma, y0) = new(mu, sigma, y0)
end


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


function collapse_but_first(A::AbstractArray)
    s = size(A)
    dims = tuple(find([d == 1 && i != 1 for (i,d) in enumerate(s)])...)
    return squeeze(A,dims)
end

vwrap(x::AbstractArray) = x
vwrap(x::Number) = [x]

@inline function add!{A}(y::A, dy::A)
  for i in eachindex(y)
    @inbounds y[i] += dy[i]
  end
end

@inline function cmul!{A}(y::A, dy::Float64)
  for i in eachindex(y)
    @inbounds y[i] *= dy
  end
end


end