using Revise, ApproxOperator, LinearAlgebra, XLSX
include("input.jl")

ndiv= 100
ndiv_p= 1
elements,nodes,nodes_p = import_fem_bar("./msh/bar_"*string(ndiv)*".msh","./msh/bar_"*string(ndiv_p)*".msh")
nᵤ = length(nodes)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
i=1/10
R = 1
h = R*i
E = 3e6
I = h^3/12
A = h
EI = E*I
EA = E*A
kGA = EA/2*5/6
P  = 1000

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

eᵇ = d₁[2]/u
eˢ = d₂[2]/v
eᵐ = d₃[2]/θ

index = 100
XLSX.openxlsx("./xlsx/curved.xlsx", mode="rw") do xf
    Sheet = xf[2]
    ind = findfirst(n->n==ndiv,index)+8
    Sheet["A"*string(ind)] = i
    Sheet["B"*string(ind)] = eᵇ
    Sheet["C"*string(ind)] = eˢ
    Sheet["D"*string(ind)] = eᵐ
end

θᶠ/θᶜ