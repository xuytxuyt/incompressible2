using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

elements,nodes,nodes_p = import_fem_tri3("./msh/square_8.msh","./msh/square_8.msh")
# elements,nodes,nodes_p= import_quad("./msh/square_quad_2.msh","./msh/square_quad_2.msh")

nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

E = 3e6
# ν=0.3
ν=0.49999999999999
Ē = E/(1-ν^2)
ν̄ = ν/(1-ν)
u(x,y) = x+y
v(x,y) = x-y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = -1.0

ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1-ν)*n₁+E/(1+ν)*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1-ν)*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)/(1-2ν)*n₁+E/(1+ν)*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1+ν)/(1-2ν)*n₂)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1+ν)*n₂)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁-E/(1+ν)*n₂)

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫qpdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:g}(),
    Operator{:Hₑ_PlaneStress}(:E=>Ē,:ν=>ν̄),
]

k = zeros(2*nₚ,2*nₚ)
kᵤ = zeros(2*nₚ,nₚ)
kₚ = zeros(nₚ,nₚ)
f = zeros(2*nₚ)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω"],k)
ops[3](elements["Ω"],k)
ops[4](elements["Ω"],elements["Ω"],kᵤ)
ops[5](elements["Ω"],kₚ)
ops[6](elements["Γᵍ"],k,f)
ops[7](elements["Γᵗ"],f)

k = [k kᵤ;kᵤ' zeros(nₚ,nₚ)]
# k = [k kᵤ;kᵤ' kₚ]
f = [f;zeros(nₚ)]

d = k\f
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

Hₑ_PlaneStress = ops[9](elements["Ω"])