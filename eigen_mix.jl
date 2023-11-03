using Revise, ApproxOperator, LinearAlgebra, Printf, TimerOutputs, SparseArrays
include("input.jl")

ndiv= 8
ndiv_p= 16
elements,nodes,nodes_p = import_fem_tri3("./msh/square_"*string(ndiv)*".msh","./msh/square_"*string(ndiv_p)*".msh")

nᵤ = length(nodes)
nₚ = length(nodes_p)

s = 2.5*10/ndiv_p*ones(nₚ)
push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
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

prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3*x)*x + (2+ν)*(y^2-D^2/4)))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫∫p∇vdxdy}(),
       Operator{:∫∫qpdxdy}(:E=>E,:ν=>ν),
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

# kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nᵤ,2*nᵤ)
kᵍ = zeros(2*nᵤ,2*nᵤ) 
kᵤₚ = zeros(2*nᵤ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
# kₚ = zeros(nᵤ,nᵤ)
f = zeros(2*nᵤ)

opsᵈ[1](elements["Ω"],kᵈ)
ops[2](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
ops[3](elements["Ωᵖ"],kₚₚ)
ops[5](elements["Γᵍ"],kᵍ,f)


# k=kᵤₚ*inv(kₚₚ)*kᵤₚ'
k=kᵤₚ/kₚₚ*kᵤₚ'

a = eigvals(k,kᵈ+kᵍ)


