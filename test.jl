using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")
elements,nodes,elms_p = import_tri3("./msh/cantilever_2.msh","./msh/cantilever_2.msh")

# nₚ = length(nodes)


# s = 2.5*12/2*ones(nₚ)

# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

# set𝝭!(elements["Ω"]) 