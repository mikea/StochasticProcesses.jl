immutable PoissonProcess
    lambda::Float64
end

type PoissonState
    t::Vector{Float64}
    n::Int64
    k::Int64
    y::Vector{Float64}
    dt::Vector{Float64}
    
    PoissonState(t, k) = new(t, length(t), k, zeros(k), diff(t))
end

function next(p::PoissonProcess, s::PoissonState, i)
    if i == 1
        return s.y
    end
    
    s.y .= s.y .+ rand(Poisson(p.lambda * s.dt[i-1]), s.k)
    s.y
end

function cumsim(p::PoissonProcess, t, k)
    s = PoissonState(t, k)
    result = zeros(s.n, s.k)
    for i in 1:s.n
        y = next(p, s, i)
        result[i,:] .= y[:]
    end
    result
end

function sim(p::PoissonProcess, t, k)
    s = PoissonState(t, k)
    for i in 1:(s.n-1)
        next(p, s, i)
    end
    next(p, s, s.n)
end

distribution(p::PoissonProcess, t) = Poisson(p.lambda * t)