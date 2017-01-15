using Base.Test
using StochasticProcesses


let p = FirstTime(BrownianMotion(), (t,x)->x>.1), t=linspace(0, 1, 10)
  sim(p, t, 3)
end