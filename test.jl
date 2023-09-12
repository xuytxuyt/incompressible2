using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

elements,nodes = import_mf_tri3("./msh/cantilever_2.msh","./msh/cantilever_2.msh")
# elements= import_mf_tri3("./msh/cantilever_quad_2.msh","./msh/cantilever_quad_2.msh")

# nₚ = length(nodes)


# s = 2.5*12/2*ones(nₚ)

# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
