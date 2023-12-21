using Revise, ApproxOperator, LinearAlgebra, XLSX
include("input.jl")

ndiv= 30
ndiv_n= 30
ndiv_v= 30
elements,nodes,nodes_n,nodes_v = import_fem_bar("./msh/bar_"*string(ndiv)*".msh","./msh/bar_"*string(ndiv_n)*".msh","./msh/bar_"*string(ndiv_v)*".msh")
nₖ = length(nodes)
nₙ = length(nodes_n)
nᵥ = length(nodes_v)

sₙ = 1.5/ndiv_n*ones(nₙ)
sᵥ = 1.5/ndiv_v*ones(nᵥ)

push!(nodes_n,:s₁=>sₙ,:s₂=>sₙ,:s₃=>sₙ)
push!(nodes_v,:s₁=>sᵥ,:s₂=>sᵥ,:s₃=>sᵥ)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωⁿ"])
set𝝭!(elements["Ωᵛ"])
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
    #    Operator{:∫κεγds}(:EI=>EI,:EA=>EA,:kGA=>kGA,:R=>R),
       Operator{:∫κMds}(:EI=>EI),
       Operator{:∫εNds}(:R=>R),
       Operator{:∫γVds}(:R=>R),
       Operator{:∫nNds}(:EA=>EA),
       Operator{:∫vVds}(:kGA=>kGA),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢθᵢds}(:α=>1e9*E),
   #     Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
   #     Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
]

kᵏᵏ = zeros(3*nₖ,3*nₖ)
kᵏⁿ = zeros(3*nₖ,nₙ)
kᵏᵛ = zeros(3*nₖ,nᵥ)
kⁿⁿ = zeros(nₙ,nₙ)
kᵛᵛ = zeros(nᵥ,nᵥ)
f = zeros(3*nₖ)
f[3*nₖ-1] += P

ops[1](elements["Ω"],kᵏᵏ)
ops[2](elements["Ω"],elements["Ωⁿ"],kᵏⁿ)
ops[3](elements["Ω"],elements["Ωᵛ"],kᵏᵛ)
ops[4](elements["Ωⁿ"],kⁿⁿ)
ops[5](elements["Ωᵛ"],kᵛᵛ)
ops[7](elements["Γᵍ"],kᵏᵏ,f)

k = [kᵏᵏ kᵏⁿ kᵏᵛ;kᵏⁿ' kⁿⁿ zeros(nₙ,nᵥ);kᵏᵛ' zeros(nᵥ,nₙ) kᵛᵛ]
f = [f;zeros(nₙ);zeros(nᵥ)]
d = k\f
d₁ = d[1:3:3*nₖ]
d₂ = d[2:3:3*nₖ]
d₃ = d[3:3:3*nₖ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)

θ = P*R^2/EI
u = P*R^3/2/EI-P*R/2/kGA-P*R/2/EA
v = π*P*R^3/4/EI+π*P*R/4/kGA+π*P*R/4/EA

eᵇ = d₁[2]/u
eˢ = d₂[2]/v
eᵐ = d₃[2]/θ
