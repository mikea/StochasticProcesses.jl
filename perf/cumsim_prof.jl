using StochasticProcesses

bm = BrownianMotion()
t = linspace(0, 1, 1000)

@time cumsim(bm, t, 1000)
Profile.clear()
@profile @time cumsim(bm, t, 1000)

f = open(perf/"cumsim_prof.prof","w")
withenv("COLUMNS"=>"10000") do
  Profile.print(f)
end
close(f)