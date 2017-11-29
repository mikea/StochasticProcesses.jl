using Base.Test
using StochasticProcesses

t=linspace(0, 1, 10)

let process = BrownianMotionWithDrift(.1, .2, 1.0)
  @test initial(process) ≈ 1.0

  # [0] time grid
  @test sim(process,    [0])      ≈ [1.0]
  @test sim(process,    [0], 3)   ≈ [1.0 1.0 1.0]
  @test cumsim(process, [0])      ≈ [1.0]
  @test cumsim(process, [0], 3)   ≈ [1.0 1.0 1.0]


  @test size(sim(process, [0, 1]))      == (1,)
  @test size(sim(process, [0, 1], 3))   == (1,3)
  @test size(cumsim(process, [0,1]))    == (2,)
  @test size(cumsim(process, [0,1], 3)) == (2,3)

  distribution(process, 0)
  distribution(process, 10)
  @test size(rand(process, t, 10)) == (10,)
end


# multivariate bm
let process = BrownianMotionWithDrift([.1, .2], eye(2,2), [1.0, 2.0])
  @test sim(process, [0]) ≈ [1.0, 2.0]
  @test sim(process, [0], 3) ≈ [1.0 1.0 1.0; 2.0 2.0 2.0]
  @test cumsim(process, [0]) ≈ [1.0 2.0]
  @test cumsim(process, [0], 3)[1,:,:] ≈ [1.0 1.0 1.0; 2.0 2.0 2.0]

  @test size(sim(process, [0, 1]))      == (2,)
  @test size(sim(process, [0, 1], 3))   == (2,3)
  @test size(cumsim(process, [0, .5, 1]))    == (3,2)
  @test size(cumsim(process, [0, .5, 1], 4)) == (3,2,4)

  distribution(process, 0)
  distribution(process, 10)

  @test size(rand(process, t, 10)) == (2, 10)
end
