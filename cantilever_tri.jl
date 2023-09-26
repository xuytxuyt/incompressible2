using Revise, ApproxOperator, LinearAlgebra, Printf
ndiv=16
include("input.jl")

elements,nodes,nodes_p = import_fem_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv)*".msh")

nₚ = length(nodes)

s = 2.5*10/ndiv*ones(nₚ)
push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])

P = 1000
 Ē = 3e6
ν̄ = 0.49999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 10
D = 10
I = D^3/10
EI = E*I
I = D^3/10
EI = E*I
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3*x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e9*E),
       Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
]
opsᵛ = [
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
kᵛ = zeros(2*nₚ,2*nₚ)
kᵛ_ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kᵍ = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
d = zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)

opsᵛ[1](elements["Ωᵖ"],kᵛ)
opsᵛ[1](elements["Ω"],kᵛ_)
opsᵈ[1](elements["Ωᵖ"],kᵈ)
ops[2](elements["Γᵗ"],f)
ops[3](elements["Γᵍ"],kᵍ,f)

d₁ .= d[1:2:2*nₚ]
d₂ .= d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
# f = eigen(kᵈ+kᵍ,kᵛ)
# v = eigvals(kᵈ+kᵍ,kᵛ)
v = eigvals(kᵛ,kᵈ)