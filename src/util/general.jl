subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))
nargs_f(f::Function)  = first(methods(f)).nargs-1 # number of arguments of f
