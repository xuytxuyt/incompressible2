
using  ApproxOperator, LinearAlgebra, Printf, XLSX

include("input.jl")

ndiv= 60
elements,nodes= import_quad_PP("./msh/cantilever_quad_"*string(ndiv)*".msh")
nᵤ = length(nodes)
nₚ = length(elements["Ωᵖ"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Ωᵍ"])
# set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])

P = 1000
Ē = 3e6
ν̄ = 0.49999999999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 48
D = 12
I = D^3/12
EI = E*I
ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
Operator{:∫∫p∇vdxdy}(),
Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
Operator{:∫vᵢtᵢds}(),
Operator{:∫vᵢgᵢds}(:α=>1e9*E),
Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    ]
kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kᵤₚ = zeros(2*nᵤ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(2*nᵤ)

ops[3](elements["Ω"],kᵤᵤ)
ops[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
ops[5](elements["Ωᵖ"],kₚₚ)
ops[7](elements["Γᵍ"],kᵤᵤ,f)
ops[6](elements["Γᵗ"],f)

k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;zeros(nₚ)]

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)

h1,l2 = ops[8](elements["Ω"])
L2 = log10(l2)
H1 = log10(h1)
h = log10(12.0/ndiv)

index = [8,16,32,60]
XLSX.openxlsx("./xlsx/eigem.xlsx", mode="rw") do xf
    Sheet = xf[3]
    ind = findfirst(n->n==ndiv,index)+1
    Sheet["B"*string(ind)] = h
    Sheet["C"*string(ind)] = L2
    Sheet["D"*string(ind)] = H1
 end

