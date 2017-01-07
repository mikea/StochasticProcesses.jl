abstract AResult{T,N} <: AbstractArray{T,N}

immutable SimResult{T,N} <: AResult{T,N}
  sol::Array{T,N}
  
  t # todo: type this
  b::Array{T, 2}
end

immutable CumsimResult{T,N} <: AResult{T,N}
  sol::Array{T,N}
  t # todo: type this
end


Base.size(r::AResult) = size(r.sol)
Base.getindex(r::AResult, i...) = getindex(r.sol, i...)
