"""Entry point for extracting linear Jacobian blocks from a JuMP model."""

module linear_jac

using JuMP
import MathOptInterface as MOI

# Return variable vector (creation order) and names
function _vars_and_names(m::Model)
    vars = all_variables(m)
    names = String[]
    sizehint!(names, length(vars))
    for v in vars
        push!(names, JuMP.name(v))
    end
    return vars, names
end

# Build a coefficient row for an AffExpr using the given variable order
function _row_from_affexpr(f::AffExpr, vars::Vector{VariableRef})
    n = length(vars)
    row = zeros(Float64, n)
    # Access coefficients directly from the AffExpr structure
    for (var, coef) in f.terms
        # Find the index of this variable in our ordered list
        idx = findfirst(==(var), vars)
        if !isnothing(idx)
            row[idx] = coef
        end
    end
    return row
end

# Push 1×n row
_push_row!(store::Vector{Vector{Float64}}, r::Vector{Float64}) = (push!(store, r); store)

"""
    extract_linear_blocks(m; obj_exprs = AffExpr[])

Extract linear Jacobian blocks from JuMP model `m`.

Returns `(Jeq, Jineq, Jobj, varnames)` where:

- `Jeq`: equality rows (for `f(x) == b`)
- `Jineq`: inequalities standardized to “≤ 0”
    * `f(x) ≤ u`           → row =  ∇f
    * `l ≤ f(x)`           → row = −∇f    (since `l − f(x) ≤ 0`)
    * `l ≤ f(x) ≤ u`       → two rows:  ∇f   and   −∇f
- `Jobj`: one row per affine objective expression in `obj_exprs`
- `varnames`: column order equals `all_variables(m)`

Assumes all constraints/expressions you care about are **linear (AffExpr)**.
"""
function extract_linear_blocks(m::Model; obj_exprs::Vector{AffExpr}=AffExpr[])
    vars, names = _vars_and_names(m)
    n = length(vars)

    eq_rows = Vector{Vector{Float64}}()
    le_rows = Vector{Vector{Float64}}()

    # Equalities: AffExpr in EqualTo
    for cref in all_constraints(m, AffExpr, MOI.EqualTo{Float64})
        con = JuMP.constraint_object(cref)
        f   = con.func
        row = _row_from_affexpr(f, vars)
        _push_row!(eq_rows, row)
    end

    for cref in all_constraints(m, AffExpr, MOI.LessThan{Float64})
        con = JuMP.constraint_object(cref)
        f   = con.func
        row = _row_from_affexpr(f, vars)
        _push_row!(le_rows, row)           # f - u ≤ 0 → ∇f
    end

    for cref in all_constraints(m, AffExpr, MOI.GreaterThan{Float64})
        con = JuMP.constraint_object(cref)
        f   = con.func
        row = _row_from_affexpr(f, vars)
        _push_row!(le_rows, -row)          # l - f ≤ 0 → -∇f
    end

    for cref in all_constraints(m, AffExpr, MOI.Interval{Float64})
        con = JuMP.constraint_object(cref)
        f   = con.func
        row = _row_from_affexpr(f, vars)
        _push_row!(le_rows,  row)          # f - u ≤ 0
        _push_row!(le_rows, -row)          # l - f ≤ 0
    end


    # Dense matrices
    Jeq   = isempty(eq_rows) ? zeros(0, n) : reduce(vcat, (permutedims(r) for r in eq_rows))
    Jineq = isempty(le_rows) ? zeros(0, n) : reduce(vcat, (permutedims(r) for r in le_rows))

    # Objective “gradients”: coefficient rows of the provided AffExprs
    Jobj_rows = Vector{Vector{Float64}}()
    for expr in obj_exprs
        push!(Jobj_rows, _row_from_affexpr(expr, vars))
    end
    Jobj = isempty(Jobj_rows) ? zeros(0, n) : reduce(vcat, (permutedims(r) for r in Jobj_rows))

    return Jeq, Jineq, Jobj, names
end

end # module
