using Base.Test
using StochasticProcesses

t=linspace(0, 1, 10)

let process = BrownianMotion()
  @test size(cumsim(process, t, 3)) == (10, 3)
  @test length(sim(process, t, 3)) == 3
  distribution(process, 0)
  distribution(process, 10)
  rand(process, t, 10)
end
