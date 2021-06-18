# unit_cell+lattice -> matrix of Hamiltonian and Current operators
using .Util.Physics
using SparseArrays, LinearAlgebra
export get_operator_gen, get_hamiltonian, get_current_op, get_operator



# returns sparse matrix generator
# functions applying to all terms:[[NOT THERE YET!!!]]
# f_mult: xloc_ann, xloc_cre -> y, and this y is *multiplied* to the term between ann and cre
# f_add: xloc_ann, xloc_cre -> y, and this y is *added* to the term between ann and cre
# f_comp: x -> y, hence this f is composed, i.e. it acts upon the value of the term.
function get_operator_gen(ltc::Lattice; 
                          nnz_compress=true, 
                          nt_compress=true,
                          idx_compress=true,
                          j_to_nnz_table=true,
                         )

    # prepare
    refresh_none(ltc)

    #    nnz, uc_nnz = est_nnz(ltc)
    NNZ_est = est_nnz(ltc)
    NT_est = get_term_count(ltc)
    dof_tot = get_dof(ltc)
    dof_per_uc_all = map(get_dof, values(ltc.uc_table))
    dof_per_uc = dof_per_uc_all[1]
    @assert all(y->y==dof_per_uc, dof_per_uc_all) "Require all uc to have identical dof for now, but getting $(dof_per_uc_all)."


    spmatgen = SparseMatrixGen{Float64, ComplexF64}(d=ltc.d, sizes=ltc.sizes, pv=ltc.pv, dof_tot=dof_tot, dof_per_uc=dof_per_uc, NT_est=NT_est, NNZ_est=NNZ_est)
    collect_f_table!(ltc, spmatgen.f_table)


    all_iloc = get_all_ind(ltc)
    filter!(x->ltc.uc_map[x...]!=:none, all_iloc);

    NT_NNZ = [1, 1]
    for iloc in all_iloc
        #        get_operator_ann_single_uc(ltc, [iloc...], I, J, V, counter, uc_n; f=identity)
        get_operator_ann_single_uc(ltc, [iloc...], spmatgen, NT_NNZ)
    end
    spmatgen.NNZ[1] = NT_NNZ[2] - 1
    spmatgen.NT[1] = NT_NNZ[1] - 1
    println("NT, NNZ = $(NT_NNZ); est $(NT_est), $(NNZ_est)")
    @assert NT_NNZ[1]<=NT_est+1 "Expecting NT_NNZ[2]<=NT_est+1 but NT_NNZ is $(NT_NNZ), NT estimation is $(NT_est)"
    @assert NT_NNZ[2]<=NNZ_est+1 "Expecting NT_NNZ[2]<=NT_est+1 but NT_NNZ is $(NT_NNZ), NT estimation is $(NNZ_est)"

    if nnz_compress
        nnz_compress!(spmatgen)
    end
    if nt_compress
        nt_compress!(spmatgen)
    end
    if idx_compress
        idx_compress!(spmatgen)
    end
    if j_to_nnz_table
        set_j_to_nnz_table!(spmatgen)
    end
    return spmatgen
    #    return sparse(I[1:counter[1]-1], J[1:counter[1]-1], V[1:counter[1]-1])
end

# J for annihilation, I for creation
#function get_operator_ann_single_uc(ltc::Lattice, iloc_ann::Array{Int,1}, I::Array{Int64,1}, J::Array{Int64,1}, V::Array{ComplexF64,1}, cnt::Array{Int64,1}, uc_n::Int64; f::Function)
function get_operator_ann_single_uc(
                                    ltc::Lattice, 
                                    iloc_ann::Array{T,1} where {T<:Integer},
                                    spmatgen::SparseMatrixGen,
                                    NT_NNZ::Array{T,1} where {T<:Integer}
                                   )
    # helper function
    get_base = x -> (spmatgen.dof_per_uc * (ijk2I(x, ltc) - 1))

    # uc_ann info
    #print("from: $iloc_ann")
    bool_ann, uc_ann_symb = get_uc_symb(ltc, iloc_ann)
    #println(", $bool_ann, $uc_ann_symb")
    if !bool_ann
        return
    end
    uc_ann = get_uc(ltc, uc_ann_symb)


    J_base = get_base(iloc_ann)

    #uc_cre_symb = uc_symb
    #uc_cre = uc_ann ### Note: now assuming current UC and target UC are the same. TODO incorporate different UC.  Will need to work on definition of t_ex. May depend on different terms

    iloc_cre = similar(iloc_ann)
    # external hopping
    for t in uc_ann.t_ex
        # info from t
        ann = t.first.atom_a # atom symbol of the annihilation operator
        cre = t.first.atom_c # atom symbol of the creation operator
        iloc_diff = t.first.uc_c # distance to the creation operator in uc count(iloc_cre - iloc_ann)
        f_symb = t.second.t0
        f_symb = Symbol(join([string(uc_ann_symb), string(f_symb)]))

        # derived info from t for annihilation operator 
        IrJridx = ((NT_NNZ[1]-1)*spmatgen.d+1):NT_NNZ[1]*spmatgen.d
        xloc_ann = view(spmatgen.Jr, IrJridx)
        #view(spmatgen.Jr_ltc, IrJridx) .= iloc_ann
#@time        xloc_ann .= get_xloc(ltc, iloc_ann, ann)
        get_xloc!(ltc, iloc_ann, ann, xloc_ann)

        # derived info from t for creation operator 
        @. iloc_cre = iloc_ann + iloc_diff
        #println(iloc_cre)
        I_base = get_base(iloc_cre)
        xloc_cre = view(spmatgen.Ir, IrJridx)
        #view(spmatgen.Ir_ltc, IrJridx) .= iloc_cre
#@time        xloc_cre .= get_xloc(ltc, iloc_cre, cre) # bool tells if the uc and operator exist TODO maybe check existence in a separate method? 
        get_xloc!(ltc, iloc_cre, cre, xloc_cre)

        # uc_cre info
        #println("uc_cre")
        bool_cre, uc_cre_symb = get_uc_symb(ltc, iloc_cre) # this will give none if outside.
        uc_cre = get_uc(ltc, uc_cre_symb)
        ##########NOTE
        #we are overwriting skipped ones
        #by not incrementing NT and NNZ counts. 
        #This is faster than taking the short circuit with if
        spmatgen.f[NT_NNZ[1]] = f_symb

        # I and J indices for new entries to be added 
        J_new = J_base .+ get_idx_range(uc_ann, ann)
        I_new = I_base .+ get_idx_range(uc_cre, cre)


        # prepare to update
        all_int_ind = get_all_ind([get_dof(uc_cre, cre),
                                   get_dof(uc_ann, ann)])
        # decide whether to increment NT and NNZ
        #bool = (uc_ann_symb!=:none) & (uc_cre_symb!=:none)
        bool = bool_cre & (uc_cre_symb!=:none) # ann is already guaranteed to exist and not none
        #println(bool)
        for int_ind in all_int_ind
            spmatgen.J[NT_NNZ[2]] = J_new[int_ind[2]] # ann
            spmatgen.I[NT_NNZ[2]] = I_new[int_ind[1]] # cre
            #            V[cnt[1]] = V_new[int_ind...]
#            print("NT_NNZ: $NT_NNZ -> ")
            NT_NNZ[2] += 1 * bool ### UPDATING NT_NNZ
        end

        spmatgen.nnz_per_term[NT_NNZ[1]] = length(all_int_ind) ##nnz

        NT_NNZ[1] += 1 * bool # NT add 1 ### UPDATING NT_NNZ
        #println("IBASE$(I_base)")
        #println(I_new)
        #println("$NT_NNZ where $(sum(spmatgen.I[1:NT_NNZ[2]].==62))")
    end
end

function dir_to_Jdir(dir::Union{Symbol, String})
    return Symbol(string(:J, dir))
end


get_hamiltonian(ltc::Lattice) = get_operator(ltc; Ops=Symbol[:H])
get_current_op(ltc::Lattice; dir::Symbol=:x, da::Array=[1.0, 0,0], db::Array=[1.0,0,0]) = get_operator(ltc; Ops=Symbol[dir_to_Jdir(dir)], da=da, db=db)
function get_operator(ltc; kwargs...)
    Hsp_gen = get_operator_gen(ltc)
    All_ops = to_sparse_matrix(Hsp_gen; kwargs...)
    return All_ops
end

f_gen_operator(op_name::Symbol; d::Integer=3, da::Array{Float64,1}=[1.0,0,0], db::Array{Float64,1}=[0,1.0,0]) = f_gen_operator(op_name, d; da=da, db=db)
function f_gen_operator(op_name::Symbol, d::Integer; 
                        da::Array{Float64, 1}=[1.0,0,0],
                        db::Array{Float64, 1}=[0.0,1,0]
                       )
    if op_name == :H
        return (r1, r2) -> 1.0
    elseif string(op_name)[1:1]=="J"
        op_dir = Symbol(string(op_name)[2:2])
        return f_gen_current_operator(op_dir, d; da=da, db=db)
    else
        println("$(op_name) not recognized.")
    end
end


# dispatch the function to generate the function multiplied on each term 
# for defining current operator
f_gen_current_operator(dir::Symbol; d::Integer=3,
                       da::Array{Float64, 1}=[1.0,0,0],
                       db::Array{Float64, 1}=[0.0,1,0]
                      ) = f_gen_current_operator(dir, d;da=da, db=db)
function f_gen_current_operator(dir::Symbol, d::Integer; 
                                da::Array{Float64, 1}=[1.0,0,0],
                                db::Array{Float64, 1}=[0.0,1,0]
                               )
    if dir==:x || dir==:X
        dir_vec = zeros(Float64, d)
        dir_vec[1] += 1
        return f_gen_current_operator(dir_vec)

    elseif dir==:y || dir==:Y
        dir_vec = zeros(Float64, d)
        dir_vec[2] += 1
        return f_gen_current_operator(dir_vec)
    elseif d>=3 && dir==:z || dir==:Z
        dir_vec = zeros(Float64, d)
        dir_vec[3] += 1
        return f_gen_current_operator(dir_vec)
    elseif dir==:a
        @assert (d==length(da)) "dimension mismatch: $(da) should have d=$(d)"
        if !(norm(da)≈1.0)
            println("da = $(da) not normalized. Normalizing ...")
            da = da / norm(da)
        end
        return f_gen_current_operator(da)
    elseif dir==:b
        @assert (d==length(db)) "dimension mismatch: $(db) should have d=$(d)"
        if !(norm(db)≈1.0)
            println("db = $(db) not normalized. Normalizing ...")
            db = db / norm(db)
        end
        return f_gen_current_operator(db)
    else
        println("Warning: current operator defined incorrectly.")
        return identity
    end
end

# generate function multiplied on each term for defining current operator
# See my notes. TODO link to notes
function f_gen_current_operator(dir_vec::Array{T,1} where {T<:Real})
    @assert norm(dir_vec)≈1.0
    f = (xloc_ann,
         xloc_cre) -> -1.0im * dot(xloc_cre - xloc_ann, dir_vec)
    return f
end


# Now checking existence depend on this step and get_uc_symb
function get_xloc(ltc::Lattice, iloc::Array{Int64,1}, atom_symb::Symbol)
    xloc = similar(ltc.pv[1])
    get_xloc!(ltc, iloc, atom_symb, xloc)
    return xloc
#    else
#        return [0.0, 0.0]
#    end
end
function get_xloc!(ltc::Lattice, iloc::Array{Int64,1}, atom_symb::Symbol, xloc::Union{Array, SubArray})
    bool, uc_symb = get_uc_symb(ltc, iloc)
    uc = get_uc(ltc, uc_symb)
    get_xloc!(ltc, iloc, xloc)
    add_xloc!(uc, atom_symb, xloc)
end


function apply_c_mult(f_table::Dict{Symbol, Function}, f_mult::Function)
    # NOT UPDATED
    new_f_table = Dict(f_table)
    for (f_symb, f) in new_f_table
        #println("create new $(f_symb) to multiply $(f) with $(f_mult)")
        g =  function(x, y, v)
                #v = f(x, y) 
                #v .*= f_mult(x, y)
                #return v
                #return f(x,y) * f_mult(x,y)
                f(x, y, v)
                v .*= f_mult(x, y)
            end
        #g = Meta.parse("(x, y) -> ((($(f))(x,y)) * (($(f_mult))(x,y)))")
        new_f_table[f_symb] = g
    end
    return new_f_table
end

function apply_c_mult!(f_table::Dict{Symbol, Function}, f_mult::Function)
    # NOT UPDATED
    for (f_symb, f) in f_table
        println("modify $(f_symb) to multiply f_mult")
        g =  function(x,y,v)
                #v = f(x, y) 
                #v .*= f_mult(x, y)
                #return v
                #return f(x,y) * f_mult(x,y)
                f(x, y, v)
                v .*= f_mult(x, y)
            end
       # g = Meta.parse("(x, y) -> ((($(f))(x,y)) * (($(f_mult))(x,y)))")
        f_table[f_symb] = g
    end
end

allsame(arr) = all(x->x==arr[1], arr)

#function eval_f(f_table)
#    f_table_evaled = Dict{Symbol, Function}()
#    for (f_symb, f) in f_table
#        f_table_evaled[f_symb] = eval(f)
#    end
#    return f_table_evaled
#end



σx, σy, σz, σ0 = sigmaMatrices(3)

