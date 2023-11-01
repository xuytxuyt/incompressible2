using Revise, ApproxOperator, LinearAlgebra, Printf
include("input.jl")

# elements,nodes,nodes_p = import_mf_tri3("./msh/square_8.msh","./msh/square_8.msh")
elements,nodes,nodes_p = import_fem_tri3("./msh/square_4.msh","./msh/square_8.msh")

nᵤ = length(nodes)
nₚ = length(nodes_p)

s = 1.5*10/8*ones(nₚ)
push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
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
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1+ν)*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁-E/(1+ν)*n₂)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)/(1-2ν)*((1-ν)*∂u∂x(x,y) + ν*∂v∂y(x,y))*n₁+E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₂)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₁+E/(1+ν)/(1-2ν)*(ν*∂u∂x(x,y) + (1-ν)*∂v∂y(x,y))*n₂)

ApproxOperator.prescribe!(elements["Γᵗ"],:g₁=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Γᵗ"],:g₂=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z)->1.0)

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

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚₚ= zeros(nₚ,nₚ)
kᵤₚ = zeros(2*nᵤ,nₚ)
f = zeros(2*nᵤ)

# ops[1](elements["Ω"],k)
# ops[2](elements["Ω"],k)
ops[3](elements["Ω"],kᵤᵤ)
# ops[4](elements["Ω"],elements["Ω"],kᵤ)
# ops[3](elements["Ωᵖ"],k)
ops[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
ops[5](elements["Ωᵖ"],kₚₚ)
ops[6](elements["Γᵍ"],kᵤᵤ,f)
# ops[7](elements["Γᵗ"],f)
ops[6](elements["Γᵗ"],kᵤᵤ,f)

k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;zeros(nₚ)]

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

Hₑ_PlaneStress = ops[9](elements["Ω"])