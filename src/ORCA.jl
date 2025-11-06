module ORCA

using JuMP
using DataFrames
using NLPModels, NLPModelsJuMP, MathOptInterface

include("linear_jac.jl")
include("correlation_strength.jl")

# Import functions from submodules
using .linear_jac
using .correlation_strength

"""
    main(model, obj_exprs)

Entry point for extracting linear Jacobian blocks from a JuMP model.
"""
function main(model, obj_exprs, num_groups=2)
    Jeq, Jineq, Jobj, varnames = linear_jac.extract_linear_blocks(model; obj_exprs)
    # println("Columns (variables): ", varnames)
    println("Jeq size   = ", size(Jeq))
    println("Jobj = ", Jobj)
    println("Jineq size = ", size(Jineq))
    println("Jobj size  = ", size(Jobj))

    results = NLPCorrStrengGenerating(Jineq, Jobj, Jeq, num_groups)
    println("NLPCorrStrengGenerating results:", results.groups)
    return results
end

export main  # expose this as the package's public entry point

end # module
