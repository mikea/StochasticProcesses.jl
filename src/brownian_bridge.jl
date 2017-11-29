struct BrownianBridge
end

function cumsim(::BrownianBridge, t, k)
    b = cumsim(BrownianMotion(), t, k)
    t1 = t[1]
    dt = t[end] - t1
    result = zeros(b)
    for i in eachindex(t)
        result[i,:] .= b[i,:] .- b[end,:] * (t[i] - t1) / dt
    end
    result
end