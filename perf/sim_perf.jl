using StochasticProcesses

bm = BrownianMotion()
t = linspace(0, 1, 10000)
sim(bm, t, 100)

println("sim(BrownianMotion) 10000x10000:")
for i in 1:5
  @time sim(bm, t, 10000)
end
