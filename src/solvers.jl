mutable struct Generator
  # specified
  t::Vector{Float64}
  n::Int64
  "number of dimensions."
  k::Int64

  # computed
  d::Tuple
  y::Matrix{Float64}
  dt::Vector{Float64}
  sdt::Vector{Float64}
  b::Matrix{Float64}
  db::Matrix{Float64}

  function Generator(t, s, k::Int)
    g      = new(t, length(t), k)
    g.d    = (s..., k)
    g.dt   = diff(t)
    g.sdt  = sqrt.(g.dt)
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

  randn!(g.db)
  cmul!(g.db, g.sdt[i-1])
  result = next!(process, g, i)
  # Important to use previous b for process computation not to have bias.
  add!(g.b, g.db)
  result
end

function cumsim{P <: AProcess}(process::P, t, k::Int=1)
  s=size(vwrap(initial(process)))
  gen = Generator(vec(t), s, k)
  result = zeros(gen.n, s..., k)
  for i in 1:gen.n
    y = next(process, gen, i)
    result[i,:,:] .= y[:,:]
  end
  CumsimResult(collapse_but_first(result), t)
end

function sim{P <: AProcess}(process::P, t, k::Int=1)
  s=size(vwrap(initial(process)))
  gen = Generator(t, s, k)
  for i in 1:(gen.n-1)
    next(process, gen, i)
  end
  SimResult(collapse_but_first(next(process, gen, gen.n)), t[end], gen.b)
end
