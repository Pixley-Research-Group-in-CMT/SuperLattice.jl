# internal data structure for lattice
export to_sparse_matrix, to_spatial_filter_vector
struct SparseMatrixGen{Ts, Tv}
    d::Int
    sizes::Array{T, 1} where {T<:Integer}
    pv::Array{Array{T,1}} where {T<:Real}
    dof_tot::Int
    dof_per_uc::Int
    NT::Array{T, 1} where {T<:Integer}
    NNZ::Array{T, 1} where {T<:Integer}
    f_table::Dict{Symbol, Function}
    I::Array{Int64, 1} 
    J::Array{Int64, 1}
    nnz_per_term::Array{Int64, 1}
    f::Array{Symbol, 1}
    Ir::Array{Ts, 1}
    Jr::Array{Ts, 1}
    #idx_to_lattice::Array{Int64, 2} # Index (as in I, J) to lattice index (UC count, such as [1,1,2])
    #Ir_ltc::Array{Ts, 1}
    #Jr_ltc::Array{Ts, 1}

    dof_compressed::Array{T, 1} where {T<:Integer}
    j_to_nnz_table::Array{Int64, 1}
    idx_decompress_table::Array{Int64, 1}
    function SparseMatrixGen{Ts, Tv}(;d::Integer=2,
                                     sizes::Array{T, 1} where {T<:Integer}=ones(d),
                                     pv::Array{Array{T,1},1} where {T<:Real}=[[Inf, Inf, Inf],[Inf, Inf, Inf],[Inf, Inf, Inf]],
                                     dof_tot::Integer,
                                     dof_per_uc::Integer,
                                     NT_est::Integer,
                                     NNZ_est::Integer
                                    ) where {Ts, Tv}
        NNZ = [0]
        NT = [0]
        f_table = Dict{Symbol, Function}()
        I = Array{Int64, 1}(undef, NNZ_est)
        J = Array{Int64, 1}(undef, NNZ_est)
        nnz_per_term = Array{Int64, 1}(undef, NT_est)
        f = Array{Symbol, 1}(undef, NT_est)
        Ir = Array{Ts, 1}(undef, NT_est * d)
        Jr = Array{Ts, 1}(undef, NT_est * d)
        #Ir_ltc = Array{Ts, 1}(undef, NT_est * d)
        #Jr_ltc = Array{Ts, 1}(undef, NT_est * d)
        #idx_to_lattice = Array{Int64, 2}(undef, d, dof_tot)

        dof_compressed = [dof_tot]
        j_to_nnz_table = Array{Int64, 1}(undef, 1)
        idx_decompress_table = Array{Int64, 1}(undef, dof_tot)
        new{Ts, Tv}(d, sizes, pv, dof_tot, dof_per_uc,
                    NT, NNZ, f_table,
                    I, J, nnz_per_term,f, Ir, Jr, 
                    #idx_to_lattice,
                    #Ir_ltc, Jr_ltc,
                    dof_compressed, j_to_nnz_table, idx_decompress_table)
    end
end

"""
Convert SparseMatrixGen to an Array of sparse matrices, where each term of the Array is defined by Ops.
Ops may include:

:H = Hamiltonian
:[JR][xyzab] = Current (J) or Position (R) operator along direction x, y, z, or da, db.
(example: 
:Jx is current operator along x; 
:Ra is position operator along ra

tbc_theta is the twist angle, range [0, 2pi]
"""
function to_sparse_matrix(spmatgen::SparseMatrixGen{Ts, Tv}; 
                          Ops::Array{Symbol, 1}=Symbol[:H],
                          da::Array{T, 1} where T<:Real =[1.0,0,0],
                          db::Array{T, 1} where T<:Real =[0.0,1,0],
                          tbc_theta::Array{T, 1} where T<:Real=[0.0,0.0,0.0],
                          save_ascii::String="", 
                          # work space arrays
                          sparse_klasttouch::Vector{Int64}=Vector{Int64}(),
                          sparse_csrrowptr::Vector{Int64}=Vector{Int64}(),
                          sparse_csrcolval::Vector{Int64}=Vector{Int64}(),
                          sparse_csrnzval::Vector{Tv}=Vector{Tv}(),
                          sparse_csccolptr::Vector{Int64}=Vector{Int64}(),
                          sparse_cscrowval::Vector{Int64}=Vector{Int64}(),
                          sparse_cscnzval::Vector{Tv}=Vector{Tv}()
                         ) where {Ts, Tv}

    da = Ts.(da)
    db = Ts.(db)
    tbc_theta = Ts.(tbc_theta)

    rets = Array{SparseMatrixCSC, 1}(undef, length(Ops))
    Ops_ri = findall(x->startswith(String(x), "R"),Ops)
    rets[Ops_ri] = to_sparse_matrix_position(spmatgen; Ops=Ops[Ops_ri], da=da, db=db)

    Ops_i = findall(x->!startswith(String(x), "R"),Ops)
    Ops = Ops[Ops_i]
    if length(Ops) == 0
        return rets
    end

    f_tables = set_f_table(spmatgen; Ops=Ops, da=da, db=db, tbc_theta=tbc_theta)


    I = @view spmatgen.I[1:spmatgen.NNZ[1]]
    J = @view spmatgen.J[1:spmatgen.NNZ[1]]
    V = Array{Tv, 1}(undef, spmatgen.NNZ[1]) # can't avoid, and should not bottleneck as we need to hold the operators in the end anyway. However, if necessary we can do each operator one after another
    #    V_current_ops = Array{Tv, 2}(undef, length(J_dir),spmatgen.NNZ[1])


    i::Int64 = 0
    ri = Array{Ts, 1}(undef, spmatgen.d)
    rj = Array{Ts, 1}(undef, spmatgen.d)
    @assert allsame(@view spmatgen.nnz_per_term[1:spmatgen.NT[1]]) "this function deals with all term with uniform square size"
    dof_per_term = Int64(round(sqrt(spmatgen.nnz_per_term[1])))
    #V_single_term = Array{Tv, 2}(undef, dof_per_term, dof_per_term)
    #V_single_term_vec_view = @view V_single_term[:] 
    
    V_blocked = reshape(V, dof_per_term, dof_per_term, spmatgen.NT[1])
    

    d = spmatgen.d
    ops_tot_dof = spmatgen.dof_compressed[1]
    coo_len = length(I)
    resize!(sparse_klasttouch, ops_tot_dof)

    resize!(sparse_csrrowptr, ops_tot_dof + 1)
    resize!(sparse_csccolptr, ops_tot_dof + 1)

    resize!(sparse_csrcolval, coo_len)
    resize!(sparse_csrnzval, coo_len)

    for (f_i, f_table) in enumerate(f_tables)
        println("f_i=$(f_i)")
        @time @sync Threads.@thread for n in 1:spmatgen.NT[1]
            ri = @view spmatgen.Ir[((n-1)*d+1):(n*d)]
            rj = @view spmatgen.Jr[((n-1)*d+1):(n*d)]
            vv = @view V_blocked[:, :, n]

            Base.invokelatest(f_table[spmatgen.f[n]], ri, rj, vv) # Work around for world age problem

        end
        if save_ascii != ""
            println(size(spmatgen.Ir))
            println(size(spmatgen.Jr))
            println(size(V_blocked))
            ascii_fn = "$(save_ascii).$(String(Ops[fi]))"
            writedlm(ascii_fn, hcat(reshape(spmatgen.Ir[1:spmatgen.d*spmatgen.NT[1]], spmatgen.d, :)',
                                    reshape(spmatgen.Jr[1:spmatgen.d*spmatgen.NT[1]], spmatgen.d, :)',
                                    real(transpose(reshape(V_blocked, dof_per_term^2, :))),
                                    imag(transpose(reshape(V_blocked, dof_per_term^2, :)))
                                   ))
        end

        println("final step timing:")
        @time rets[f_i] = SparseArrays.sparse!(
                                         I, J, V, 
                                         ops_tot_dof, ops_tot_dof, +, sparse_klasttouch,
                                         sparse_csrrowptr, sparse_csrcolval, sparse_csrnzval,
                                         sparse_csccolptr, sparse_cscrowval, sparse_cscnzval, # V is no longer needed afterthis, so reused as output
                                        )
        dropzeros!(rets[f_i])
    end
    return rets
end

# check all functions are defined
function validate(spmatgen::SparseMatrixGen) 
    println("Validating...")
    req_keys = Set(unique(spmatgen.f[1:spmatgen.NT[1]]))
    all_keys = Set(keys(spmatgen.f_table))

    if req_keys <= all_keys
        println("all function loaded")
        return true
    else
        println("missing: $(setdiff(req_keys, all_keys)). \n fulfilled: $(all_keys)")
        return false
    end
end


function set_f_table(
                     spmatgen::SparseMatrixGen{Ts, Tv}; 
                     Ops::Array{Symbol, 1}=Symbol[:H],
                     da::Array{T, 1} where T<:Number=[1.0,0,0],
                     db::Array{T, 1} where T<:Number=[0.0,1,0],
                     tbc_theta::Array{T, 1} where T<:Number=[0.0,0.0,0.0]
                    ) where {Ts, Tv}
    da = Ts.(da)
    db = Ts.(db)
    tbc_theta = Ts.(tbc_theta)

    f_mults = map(op -> f_gen_operator(op; d=spmatgen.d, da=da, db=db), Ops)
    f_tables = map(f_mult -> apply_c_mult(spmatgen.f_table, f_mult), f_mults)
    if !all(x->x==0, tbc_theta)
        f_mult_tbc = gen_tbc(tbc_theta, spmatgen.sizes, spmatgen.pv)
        map!(f_table -> apply_c_mult(f_table, f_mult_tbc), f_tables, f_tables)
    end
    return f_tables
    #    f_tables = map(eval_f, f_tables)
end

function nnz_compress!(spmatgen::SparseMatrixGen)
    resize!(spmatgen.I, spmatgen.NNZ[1])
    sizehint!(spmatgen.I, spmatgen.NNZ[1])
    resize!(spmatgen.J, spmatgen.NNZ[1])
    sizehint!(spmatgen.J, spmatgen.NNZ[1])
end
function nt_compress!(spmatgen::SparseMatrixGen)
    resize!(spmatgen.Ir, spmatgen.NT[1]*spmatgen.d)
    sizehint!(spmatgen.Ir, spmatgen.NT[1]*spmatgen.d)
    resize!(spmatgen.Jr, spmatgen.NT[1]*spmatgen.d)
    sizehint!(spmatgen.Jr, spmatgen.NT[1]*spmatgen.d)
    resize!(spmatgen.nnz_per_term, spmatgen.NT[1])
    sizehint!(spmatgen.nnz_per_term, spmatgen.NT[1])
    resize!(spmatgen.f, spmatgen.NT[1])
    sizehint!(spmatgen.f, spmatgen.NT[1])
end

function idx_compress!(spmatgen::SparseMatrixGen)
    idx_all = unique([spmatgen.I; spmatgen.J])
    idx_compress_table = Dict{Int64,Int64}(idx_all .=> 1:length(idx_all))

    for i in 1:spmatgen.NNZ[1]
        spmatgen.I[i] = idx_compress_table[spmatgen.I[i]]
        spmatgen.J[i] = idx_compress_table[spmatgen.J[i]]
    end
    spmatgen.dof_compressed[1] = length(idx_all)
    set_j_to_nnz_table!(spmatgen)
    resize!(spmatgen.idx_decompress_table, length(idx_all))
    sizehint!(spmatgen.idx_decompress_table, length(idx_all))
    for idx_pair in idx_compress_table
        spmatgen.idx_decompress_table[last(idx_pair)] = first(idx_pair)
    end
end

function set_j_to_nnz_table!(spmatgen::SparseMatrixGen)
    resize!(spmatgen.j_to_nnz_table, spmatgen.dof_compressed[1])
    sizehint!(spmatgen.j_to_nnz_table, spmatgen.dof_compressed[1])
    for i in Base.OneTo(spmatgen.NNZ[1])
        spmatgen.j_to_nnz_table[spmatgen.J[i]] = i
    end
end



function idx_compress2!(spmatgen::SparseMatrixGen)
    tmp_I = copy(spmatgen.I)
    dof_comp = 1
    for i in eachindex(tmp_I)
        spmatgen.I[i] = dof_comp
        first_occur = findfirst(x->x==tmp_I[i], tmp_I)
        spmatgen.I[i] = min(dof_comp, spmatgen.I[first_occur])
        dof_comp += (dof_comp==spmatgen.I[i])
    end

    tmp_J = tmp_I;
    tmp_J .= spmatgen.J
    for i in eachindex(tmp_J)
        spmatgen.J[i] = dof_comp
        first_occur = findfirst(x->x==tmp_J[i], tmp_J)
        spmatgen.J[i] = min(dof_comp, spmatgen.J[first_occur])
        dof_comp += (dof_comp==spmatgen.J[i])
    end

end

"""
    to_lattice_filter_vector(spmatgen, ltc, invariant_f, invariant_val; epsilon)

Create a filter based on latex index constraint.

Inputs:
    spmatgen
    ltc
    constr_f
    constr_val

Output:
    Int64 Array with same dimension as Hamiltonian such that index of the unit cell 
    of the point satisfies `abs(constr_f(x) - constr_val) < epsilon`.
"""
function to_lattice_filter_vector(
                                  spmatgen::SparseMatrixGen{Ts, Tv},
                                  ltc::Lattice,
                                  constr_f::Function,
                                  constr_val::Number;
                                  epsilon::Number=1e-4
                                 ) where {Ts, Tv}
    function idx_to_psi(I)
        expanded_idx = I2ijk(I, ltc)
        predicate = Int64((abs(constr_f(expanded_idx) - constr_val) < epsilon))
        return predicate
    end
    return map(idx_to_psi, spmatgen.idx_decompress_table)
end

"""
    to_spatial_filter_vector(spmatgen, invariant_f, invariant_val)

Create a filter based on spatial constraint.

Inputs:
    spmatgen
    constr_f
    constr_val

Output:
    Int64 Array with same dimension as Hamiltonian such that points x satisfies
    `abs(constr_f(x) - constr_val) < epsilon`
"""
function to_spatial_filter_vector(
                                  spmatgen::SparseMatrixGen{Ts, Tv},
                                  constr_f::Function,
                                  constr_val::Number;
                                  epsilon::Number=1e-4
                                 ) where {Ts, Tv}
    @assert allsame(@view spmatgen.nnz_per_term[1:spmatgen.NT[1]]) "this function deals with all term with uniform square size"
    nnz_per_term = spmatgen.nnz_per_term[1]
    d = spmatgen.d

    try 
        tmp = constr_f(zeros(Tv, d))
        typeof(tmp) <: Number
    catch
        throw(ArgumentError("`constr_f` should take array of $(d) and output a number."))
    end

    
    ψint = zeros(Int64, spmatgen.dof_compressed[1])
    J_unique = unique(spmatgen.J[1:spmatgen.NNZ[1]])
    @assert (maximum(J_unique) == length(J_unique) == spmatgen.dof_compressed[1]) "Maximum of J should be the same as unique number of J, and same as dof_compressed. Otherwise it is not properly compressed."
    for J in J_unique
        J_NNZ = spmatgen.j_to_nnz_table[J] # find corresponding index among NNZ
        nt = Int64(ceil(J_NNZ / nnz_per_term)) # find corresponding index among NT
        r = view(spmatgen.Jr, ((nt-1)*d+1):(nt*d))
        ψint[J] |= (abs(constr_f(r) - constr_val) <= epsilon)
    end

    return ψint
end

# TODO doc
function to_sparse_matrix_position(spmatgen::SparseMatrixGen{Ts, Tv}; 
                                   Ops::Array{Symbol, 1}=Symbol[:Rx],
                                   da::Array{Float64, 1}=[1.0,0,0],
                                   db::Array{Float64, 1}=[0.0,1,0]
                                  ) where {Ts, Tv}
    @assert allsame(@view spmatgen.nnz_per_term[1:spmatgen.NT[1]]) "this function deals with all term with uniform square size"
    nnz_per_term = spmatgen.nnz_per_term[1]
    if length(spmatgen.j_to_nnz_table)==1
        println("setting j_to_nnz_table")
        set_j_to_nnz_table!(spmatgen)
    end

    d = spmatgen.d
    rets = Array{SparseMatrixCSC, 1}(undef, length(Ops))
    for i in eachindex(Ops)
        # parse direction
        dr = parse_dir(Ops[i]; d=spmatgen.d, da=da, db=db)
        
        # all based on J, the annihilation side, as the creation side may
        # have un-wrapped location.
        J_unique = unique(spmatgen.J[1:spmatgen.NNZ[1]])
        r_all = similar(J_unique, Ts)
        for (Ji, J) in enumerate(J_unique)
            J_NNZ = spmatgen.j_to_nnz_table[J]
            nt = Int64(ceil(J_NNZ / nnz_per_term))
            r = view(spmatgen.Jr, ((nt-1)*d+1):(nt*d))
            r_all[Ji] = dot(dr, r)
        end
        rets[i] = sparse(J_unique, J_unique, r_all)
    end
    return rets
end

# TODO check if this is better than `to_sparse_matrix_position`.
function to_sparse_matrix_position2(spmatgen::SparseMatrixGen{Ts, Tv}; 
                                    Ops::Array{Symbol, 1}=Symbol[:Rx],
                                    da::Array{Float64, 1}=[1.0,0,0],
                                    db::Array{Float64, 1}=[0.0,1,0]
                                   ) where {Ts, Tv}
    @assert allsame(@view spmatgen.nnz_per_term[1:spmatgen.NT[1]]) "this function deals with all term with uniform square size"
    d = spmatgen.d
    rets = Array{SparseMatrixCSC, 1}(undef, length(Ops))
    for i in eachindex(Ops)
        # parse direction
        dr = parse_dir(Ops[i]; d=spmatgen.d, da=da, db=db)
        
        # all based on J, the annihilation side, as the creation side may
        # have un-wrapped location.
        J_unique = unique(spmatgen.J[1:spmatgen.NNZ[1]])
        r_all = similar(J_unique, Ts)
        for (Ji, J) in enumerate(J_unique)
            J_first = findfirst(x->x==J, spmatgen.J)
            n = Int64(ceil(J_first / 36))
            r = view(spmatgen.Jr, ((n-1)*d+1):(n*d))
            r_all[Ji] = dot(dr, r)
        end
        rets[i] = sparse(J_unique, J_unique, r_all)
    end
    return rets
end


# TODO doc
"""
Interpret the direction. 
"""
function parse_dir(Op::Symbol;
                   d::Int64=3, 
                   da::Array{Float64, 1}=[1.0,0,0],
                   db::Array{Float64, 1}=[0.0,1,0]
                  )

    dir_str = String(Op)[end:end]

    dr = similar(da, d)
    if dir_str == "x"
        dr[1] = 1
    elseif dir_str == "y"
        dr[2] = 1
    elseif dir_str == "z"
        dr[3] = 1
    elseif dir_str == "a"
        dr .= da
    elseif dir_str == "b"
        dr .= db
    end
    return dr
end
