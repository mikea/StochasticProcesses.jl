using StochasticProcesses

bm = BrownianMotion()
t = linspace(0, 1, 10000)
cumsim(bm, t, 100)[end,:]

println("cumsim(BrownianMotion) 10000x10000:")
for i in 1:5
  @time cumsim(bm, t, 10000)[end,:]
end
