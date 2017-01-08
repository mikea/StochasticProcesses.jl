immutable FirstTime{Process, Pred}
    process::Process
    pred::Pred
end

function sim{Process, Pred}(ft::FirstTime{Process, Pred}, t, k)
  s=size(vwrap(initial(ft.process)))
  gen = Generator(t, s, k)

  times=Vector{Float64}()
  mask=falses(k)

  for i in 1:gen.n
    next(ft.process, gen, i)
    t = gen.t[i]

    # todo: we can reduce k after predicate hits.
    for j in 1:gen.k
      if !mask[j] && ft.pred(t, gen.y[j])
        mask[j] = true
        push!(times, t)
      end
    end
  end

  times
end