using ApproxOperator
include("input_gmsh.jl")

ndiv = 3
elements,nodes,elms = import_curved_beam("./msh/bar_"*string(ndiv)*".msh");

nₚ = length(nodes)
s = 1.5*π/2/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ωᵥ"])
set𝝭!(elements["Ωₙₙ"])
set𝝭!(elements["Ωₘₘ"])
set𝝭!(elements["Ωₙᵥ"])
set𝝭!(elements["Ωₘᵥ"])
set𝝭!(elements["Γᵥ"])
set𝝭!(elements["Γₙ"])
set𝝭!(elements["Γₘ"])
set𝝭!(elements["Γᵛᵥ"])
set𝝭!(elements["Γᵛₘ"])
set𝝭!(elements["Γᶿ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵐ"])
