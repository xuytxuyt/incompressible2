using Revise, ApproxOperator, LinearAlgebra, Printf, TimerOutputs, SparseArrays
ndiv=64
include("input.jl")

# elements,nodes,nodes_p = import_fem_tri3_GI1("./msh/square_"*string(ndiv)*".msh","./msh/square_"*string(ndiv)*".msh")
# elements,nodes,nodes_p = import_quad_GI1("./msh/square_quad_"*string(ndiv)*".msh","./msh/square_quad_"*string(ndiv)*".msh")
elements,nodes,nodes_p = import_quad8_GI1("./msh/square_quad8_"*string(ndiv)*".msh","./msh/square_quad8_"*string(ndiv)*".msh")

const to = TimerOutput()

nₚ = length(nodes)

@timeit to "shape function" begin
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵛ"])
set∇𝝭!(elements["Ωᵛ"])
set𝝭!(elements["Γᵍ"])
end
P = 1000
 Ē = 3e6
ν̄ = 0.49999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 12
D = 12
I = D^3/12
EI = E*I
I = D^3/12
EI = E*I

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

@timeit to "assembly matrix" begin
kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kᵍ = zeros(2*nₚ,2*nₚ) 
f = zeros(2*nₚ)


opsᵛ[1](elements["Ωᵛ"],kᵛ)
# opsᵛ[1](elements["Ω"],kᵛ)
opsᵈ[1](elements["Ω"],kᵈ)
ops[3](elements["Γᵍ"],kᵍ,f)
end

@timeit to "eigen" begin
v = eigvals(kᵛ,kᵈ+kᵍ)
end

show(to)