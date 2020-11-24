###
#  Peierls substitution
function gen_peierls(B::Array{T, 1} where {T<:Real}, d::Integer)
#    @assert (length(B)==3) "Only Hopping terms have finite peierls substitution"
    if d == 3
        f = begin
            @inline function (ra,
                      rc)
                dr = rc - ra
                l = sqrt(dot(dr, dr)) # length
                μ = dr / l # direction

                θ = 0
                # (By z0 μx + Bz x0 μy + Bx y0 μz) * l
                θ += sum(B[[2,3,1]] .* ra[[3,1,2]] .* μ) * l 

                # (By μz μx + Bz μx μy + Bx μy μz) * l^2/2
                θ += sum(B[[2,3,1]] .* μ[[3,1,2]] .* μ) * (l^2/2)

                return exp.(1.0im * θ)
            end
        end
        return f
    elseif d == 2
        f = begin
            function (ra,
                      rc)
                dr = rc - ra
                l = sqrt(dot(dr, dr)) # length
                μ = dr / l # direction

                θ = 0
                # (Bz x0 μy) * l
                θ += (B[3] * ra[1] * μ[2]) * l 

                # (Bz μx μy) * l^2/2
                θ += (B[3] * μ[1] * μ[2]) * (l^2/2)

                return exp.(1.0im * θ)
            end
        end
        return f
    end
end
gen_peierls(B::Array{T,1} where {T<:Real}; d::Integer=3) = gen_peierls(B, d)
gen_peierls(;B::Array{T,1} where {T<:Real}=[0,0,0.0], d::Integer=3) = gen_peierls(B, d)

###
# Twisted Boundary Condition
"""
generate the multiplicative factor that takes care of twisted boundary condition (TBC).
pv: primitive vector [v1, v2, v3]
"""
function gen_tbc(
                 theta::Array{T,1} where {T<:Real},
                 L::Union{Array{T,1}, T} where {T<:Integer},
                 pv::Array{Array{T, 1},1} where {T<:Real}
                 )
    @assert (!any(isinf, vcat(pv...))) "some pv is inf. "
    pv_mat = hcat(pv...)
    pv_mat_inv = inv(pv_mat)
    println("pvmatinv:$pv_mat_inv")
    println("pvmat:$pv_mat")
    θ = theta ./ L
    println("theta / L $(θ)")

    f = begin
        @inline function (ra, rc)
            exponent::Float64 = 0.0
            for i in eachindex(rc)
                for j in eachindex(theta)
                    exponent += θ[i] * pv_mat_inv[i,j] * (rc[j]-ra[j])
                end
            end

            return exp(1.0im * 2pi * exponent)
        end
    end

end

