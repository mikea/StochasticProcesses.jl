include("bm_test.jl")
include("bm_drift_test.jl")

using Base.Test
using StochasticProcesses

t=linspace(0, 1, 10)

let process = GeometricBrownianMotion(1., .1, 100)
  @test size(cumsim(process, t, 3)) == (10, 3)
  @test length(sim(process, t, 3)) == 3
  distribution(process, 0)
  distribution(process, 10);
  rand(process, t, 10)
end

let process = ItoIntegral(x->x)
  @test size(cumsim(process, t, 3)) == (10, 3)
  @test length(sim(process, t, 3)) == 3
  rand(process, t, 10)
end

let process = SDE((x,y)->x, (x,y)->y, 0.)
  @test size(cumsim(process, t, 3)) == (10, 3)
  @test length(sim(process, t, 3)) == 3
  rand(process, t, 10)
end


let process = CompositeProcess(BrownianMotion(), (t, y) -> y.^2)
  @test size(cumsim(process, t, 3)) == (10, 3)
  @test length(sim(process, t, 3)) == 3
  rand(process, t, 10)
end