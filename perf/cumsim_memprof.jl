using StochasticProcesses

bm = BrownianMotion()
t = linspace(0, 1, 10000)
cumsim(bm, t, 1000)

Profile.clear()

println("cumsim(BrownianMotion) 1000x1000:")
for i in 1:5
  @time cumsim(bm, t, 1000)
end