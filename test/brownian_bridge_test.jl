using Base.Test
using StochasticProcesses


let p = BrownianBridge(), t=linspace(0, 1, 10)
  cumsim(p, t, 3)
end