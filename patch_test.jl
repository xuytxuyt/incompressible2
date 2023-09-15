using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

elements,nodes,nodes_p = import_fem_tri3("./msh/cantilever_2.msh","./msh/cantilever_2.msh")
# elements,nodes,nodes_p= import_quad("./msh/cantilever_quad_2.msh","./msh/cantilever_quad_2.msh")

nₚ = length(nodes)
nᵖ = length(nodes_p)

s = 2.5*10/2*ones(nₚ)
sₚ = 2.5*10/2*ones(nᵖ)

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵍ"])

E = 3e6
ν=0.3
u(x,y) = x+y
v(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0

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
ApproxOperator.prescribe!(elements["Ωᵖ"],:u=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ωᵖ"],:v=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Ωᵖ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
ApproxOperator.prescribe!(elements["Ωᵖ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
ApproxOperator.prescribe!(elements["Ωᵖ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
ApproxOperator.prescribe!(elements["Ωᵖ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:g}(),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
# k = zeros(2*nᵖ,2*nᵖ)
# f = zeros(2*nᵖ)

ops[1].(elements["Ω"];k=k)
# ops[1].(elements["Ωᵖ"];k=k)
ops[2].(elements["Γᵍ"];k=k,f=f)

d = k\f
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
# d₁ = d[1:2:2*nᵖ]
# d₂ = d[2:2:2*nᵖ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
# push!(nodes_p,:d₁=>d₁,:d₂=>d₂)
Hₑ_PlaneStress = ops[4](elements["Ω"])