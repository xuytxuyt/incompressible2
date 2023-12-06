using Revise, ApproxOperator, LinearAlgebra
include("input.jl")

    ndiv= 1
    ndiv_p= 1
    elements,nodes,nodes_p = import_fem_bar("./msh/bar_"*string(ndiv)*".msh","./msh/bar_"*string(ndiv_p)*".msh")
    nᵤ = length(nodes)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

    E = 3e6
    ν = 0.3
    # ν=0.49999999999999
    G  = E/2
    I  = 1
    A  = 1    
    EI = E*I
    EA = E*A
    GA = G*A
    R  = 1
    P  = 1 
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    ops = [
           Operator{:∫κεγds}(:EI=>EI,:EA=>EA,:GA=>GA,:R=>R),
           Operator{:∫vᵢtᵢds}(),
           Operator{:∫vᵢgᵢds}(:α=>1e9*E),
           Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
           Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
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

