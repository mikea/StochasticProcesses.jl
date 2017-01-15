using Base.Test
using StochasticProcesses


let p = PoissonProcess(5), t=linspace(0, 1, 10)
  sim(p, t, 3)
  cumsim(p, t, 3)
  distribution(p, 10)
end