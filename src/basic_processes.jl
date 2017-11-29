"Brownian motion process."
struct BrownianMotion{Y} <: AItoProcess
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

struct BrownianMotionWithDrift{Y,S} <: AItoProcess
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

solution(p::BrownianMotionWithDrift, t, b) =
  p.y0 + t * p.mu + p.sigma * b

# Geometric Brownian Motion

mutable struct GeometricBrownianMotion <: AItoProcess
    mu::Float64
    sigma::Float64
    y0::Float64
end

convert(::Type{ItoProcess}, bm::GeometricBrownianMotion) = 
    ItoProcess((t,dt,b,db,y)-> (bm.mu * dt) * y + bm.sigma * (y .* db), bm.y0)

distribution(bm::GeometricBrownianMotion, t) = 
    t == 0 ? Constant(bm.y0) :
    LogNormal(log(bm.y0) + t * (bm.mu - bm.sigma^2 / 2), bm.sigma * sqrt(t))

solution(p::GeometricBrownianMotion, t, b) = 
    p.y0 * exp(t * (p.mu - p.sigma^2 / 2) + p.sigma * b)


"Ito integral of f(t, B)"
struct ItoIntegral{F} <: AItoProcess
  f::F
end

convert{F}(::Type{ItoProcess}, ii::ItoIntegral{F}) = 
    ItoProcess((t,dt,b,db,y)->ii.f(t, b).*db, 0.)

"Process corresponding to solutions of stochastic differential equation
 dy = f(t, y) dt + g(t,y) db."
struct SDE{F,G} <: AItoProcess
  f::F
  g::G
  y0::Float64
end

convert{F,G}(::Type{ItoProcess}, sde::SDE{F,G}) = 
    ItoProcess((t,dt,b,db,y)->sde.f(t, y).*dt + sde.g(t, y).*db, sde.y0)
