# ORCA
Multi-objective dimensionality reduction algorithm: Objective reduction community algorithm (ORCA).

ORCA is a Julia package designed for extracting linear Jacobian blocks and analyzing objective correlation structures in multi-objective optimization models built with JuMP.

It provides a modular interface to compute correlations between multiple objectives, making it particularly useful for dimensionality reduction and decision-making in many-objective optimization problems.

# 🧩Prerequisites

Before using ORCA, ensure you have the following:

Julia
 installed (version 1.9+ recommended)

Internet access to install package dependencies from the Julia General Registry

# ⚙️ Installation

Open a terminal (or Julia REPL) and enter:

```
julia> using Pkg
julia> Pkg.add(url="https://github.com/beam4869/ORCA.jl")
```

This will download and install the package and its dependencies.

# ✅ Verifying Installation
Once installed, you can verify everything is working by running the built-in tests:
```
julia> using Pkg
julia> Pkg.test("ORCA")
```

If everything is configured correctly, you should see an output similar to:
```
Precompiling project...
  82 dependencies successfully precompiled in 55 seconds. 110 already precompiled.
     Testing Running tests...
     Testing ORCA tests passed

```

If the tests pass — congratulations! ORCA is installed and ready to use.

# 🚀 Usage
After successful installation and testing, you can start using ORCA as follows:

```
using ORCA
using JuMP

# Example: build a JuMP model
m = Model()
@variable(m, x >= 0)
@variable(m, y >= 0)
@variable(m, z >= 0)
@constraint(m, x + y -z <= 1)
@constraint(m, 5x + y + 2z== 3)
@constraint(m, -x + y +2z <= 0)
@expression(m, obj1, x + 2y +z)
@expression(m, obj2, -3x - y -z)
@expression(m, obj3, x - y +z)

# Call ORCA's main analysis function
res = ORCA.main(m, [obj1, obj2, obj3])
println("ORCA results:", res.groups)
```

## Function signature

```
ORCA.main(model, [objective_1, objective_2, ..., objective_n], number_of_groupings)
```

`model` — a JuMP model containing your variables, constraints, and objectives

`[objective_1, objective_2, ..., objective_n]` — a list of objective expressions defined in JuMP

The number of objectives must be ≥ 2

`number_of_groupings` is default set to be 2.

# 📘 Defining Objectives in JuMP
Objective expressions in ORCA are standard JuMP expressions.
You can read more about how to define and manipulate them here: 🔗 [JuMP Expressions Documentation](https://jump.dev/JuMP.jl/stable/manual/expressions/)
