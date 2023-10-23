using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

elements,nodes,nodes_p = import_fem_tri3_direct("./msh/square_2.msh","./msh/square_2.msh")


nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

E = 3e6
ν=0.3
u(x,y) = x+y
v(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0

ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1-ν)*n₁+E/(1+ν)*n₂)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1-ν)*n₂)

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:g₂}(),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
]

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
for ap in elements["Γᵍ"]
    x, = ap.𝓒
    x.d₁ = u(x.x,x.y)
    x.d₂ = v(x.x,x.y)
end

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](elements["Ω"],k)
ops[3](elements["Γᵗ"],f)
ops[4].(elements["Γᵍ"],k=k,f=f,dof=:d₁)
ops[4].(elements["Γᵍ"],k=k,f=f,dof=:d₂)

d = k\f
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

Hₑ_PlaneStress = ops[5](elements["Ω"])