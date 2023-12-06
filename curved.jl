using Revise, ApproxOperator, LinearAlgebra
include("input.jl")

    ndiv= 2
    ndiv_p= 2
    elements,nodes,nodes_p = import_fem_bar("./msh/bar_"*string(ndiv)*".msh","./msh/bar_"*string(ndiv_p)*".msh")
    nᵤ = length(nodes)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

    E = 3e6
    EI = 3e6
    EA  = 3e6
    kGA  = EI/2*5/6
    R  = 1
    P  = 1 
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₃=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    ops = [
           Operator{:∫κεγds}(:EI=>EI,:EA=>EA,:kGA=>kGA,:R=>R),
           Operator{:∫vᵢtᵢds}(),
           Operator{:∫vᵢgᵢdΓ}(:α=>1e9*E),
       #     Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
       #     Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
    ]

    k = zeros(3*nᵤ,3*nᵤ)
    f = zeros(3*nᵤ)
    f[3*nᵤ-1] += 1.0

    ops[1](elements["Ω"],k)
    ops[3](elements["Γᵍ"],k,f)

    d = k\f
    d₁ = d[1:3:3*nᵤ]
    d₂ = d[2:3:3*nᵤ]
    d₃ = d[3:3:3*nᵤ]

    push!(nodes,:d₁=>d₁,:d₂=>d₂)

    # h1,l2 = ops[5](elements["Ω"])
    # L2 = log10(l2)
    # H1 = log10(h1)

