using  ApproxOperator, LinearAlgebra, Printf

include("input_mix.jl")

elements,nodes= import_beam("./msh/bar_test.msh")

nₚ = length(nodes)
s = 1.5*π/2/2*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω_NN"])
set𝝭!(elements["Ω_Nv"])
set𝝭!(elements["Ω_Nu"])
set𝝭!(elements["Ω_Nb"])
set𝝭!(elements["Ω_MM"])
set𝝭!(elements["Ω_Mv"])
# set𝝭!(elements["Ω_Mu"])
set𝝭!(elements["Ω_Mb"])
# set𝝭!(elements["Ω_Mθ"])
set𝝭!(elements["Ω_vN"])
set𝝭!(elements["Ω_vM"])
# set𝝭!(elements["Ω_vu"])
# set𝝭!(elements["Ω_vθ"])
# set𝝭!(elements["Ω_vt"])
# set𝝭!(elements["Ω_vm"])
set𝝭!(elements["Ω_vb"])
set𝝭!(elements["Γᵍ"])
