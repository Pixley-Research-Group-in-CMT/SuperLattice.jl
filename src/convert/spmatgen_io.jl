# save/load SparseMatrixGen
# save the generator, with anonymous function all removed, to JLD2 format
using DelimitedFiles
using JLD2, FileIO
export loadSMG, saveSMG



# remove anonymous function definition and save
function saveSMG(spmatgen::SparseMatrixGen, fn::String)
    println("removing anonymous functions for saving...")
    for k in keys(spmatgen.f_table)
        println("removing $(k)")
        delete!(spmatgen.f_table, k)
    end
    
    println("saving to $(fn)...")
    @save fn spmatgen
    println("Done saving.")
end

# load and validate functions
function loadSMG(fn::String; internal=false)
    println("reading sparse matrix generator from $(fn)...")
    @load fn spmatgen
    println("sparse matrix generator loaded.")
    if !internal
        println("call `pull_anonymous_funcs!(spmatgen, ltc)` to import functions from ltc into spmatgen")
    end
    return spmatgen
end

function loadSMG(fn::String, ltc::Lattice)
    spmatgen = loadSMG(fn; internal=true)
    pull_anonymous_funcs!(spmatgen, ltc) 
    return spmatgen
end

# pull function definition from ltc
# This is equivalent to collect_f_table but added a check.
function pull_anonymous_funcs!(spmatgen::SparseMatrixGen, ltc::Lattice)
    collect_f_table!(ltc, spmatgen.f_table)
    @assert validate(spmatgen::SparseMatrixGen) "missing function defs"
end 


