using Revise, ApproxOperator, LinearAlgebra
include("input.jl")

ndiv= 50
ndiv_p= 1
elements,nodes,nodes_p = import_fem_bar("./msh/bar_"*string(ndiv)*".msh","./msh/bar_"*string(ndiv_p)*".msh")
nᵤ = length(nodes)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

E = 3e6
EI = 3e6
EA  = 3e6
kGA  = EI/2*5/6
R  = 1
P  = 100
ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₃=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
       Operator{:∫κεγds}(:EI=>EI,:EA=>EA,:kGA=>kGA,:R=>R),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢθᵢds}(:α=>1e9*E),
   #     Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
   #     Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
]

k = zeros(3*nᵤ,3*nᵤ)
f = zeros(3*nᵤ)
f[3*nᵤ-1] += P

ops[1](elements["Ω"],k)
ops[3](elements["Γᵍ"],k,f)

d = k\f
d₁ = d[1:3:3*nᵤ]
d₂ = d[2:3:3*nᵤ]
d₃ = d[3:3:3*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)

θ = P*R^2/EI
u = P*R^3/2/EI-P*R/2/kGA-P*R/2/EA
v = π*P*R^3/4/EI+π*P*R/4/kGA+π*P*R/4/EA

eᵇ = (d₁[2]^2/u^2)^0.5
eˢ = (d₂[2]^2/v^2)^0.5
eᵐ = (d₃[2]^2/θ^2)^0.5
h = log10(1/ndiv)