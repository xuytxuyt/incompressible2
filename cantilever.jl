using Revise, ApproxOperator, LinearAlgebra, Printf, TimerOutputs, XLSX
include("input.jl")

# for i in 40:50
    ndiv= 50
    ndiv_p= 50
    # elements,nodes,nodes_p = import_fem_tri3("./msh/square_"*string(ndiv)*".msh","./msh/square_"*string(ndiv_p)*".msh")
    elements,nodes,nodes_p = import_quad("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
    # elements,nodes,nodes_p= import_quad("./msh/square_quad_"*string(ndiv)*".msh","./msh/square_quad_"*string(ndiv_p)*".msh")
    nᵤ = length(nodes)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

    P = 1000
    Ē = 3e6
    # ν̄ = 0.49999999999999
    ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 48
    D = 12
    I = D^3/12
    EI = E*I

    ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
    ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
    ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
    ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    ops = [
        Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
        Operator{:∫∫vᵢbᵢdxdy}(),
        Operator{:∫vᵢtᵢds}(),
        Operator{:∫vᵢgᵢds}(:α=>1e9*E),
        Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
        Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
    ]

    k = zeros(2*nᵤ,2*nᵤ)
    f = zeros(2*nᵤ)

    ops[1](elements["Ω"],k)
    ops[3](elements["Γᵗ"],f)
    ops[4](elements["Γᵍ"],k,f)

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]

    push!(nodes,:d₁=>d₁,:d₂=>d₂)

    h1,l2 = ops[5](elements["Ω"])
    L2 = log10(l2)
    H1 = log10(h1)
    # h = log10(10.0/ndiv)

#     index = 40:50
#     XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
#         Sheet = xf[2]
#         ind = findfirst(n->n==ndiv,index)+1
#         Sheet["F"*string(ind)] = h
#         Sheet["G"*string(ind)] = L2
#         Sheet["H"*string(ind)] = H1
#     end
# end
