using Revise, ApproxOperator, LinearAlgebra, Printf, TimerOutputs, XLSX
include("input.jl")

for i in 10:40
    ndiv= 20
    ndiv_p= i
    elements,nodes,nodes_p = import_fem_tri3("./msh/square_"*string(ndiv)*".msh","./msh/square_"*string(ndiv_p)*".msh")
    
    nᵤ = length(nodes)
    nₚ = length(nodes_p)

    s = 2.5*10/ndiv_p*ones(nₚ)
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

    P = 1000
    Ē = 3e6
    ν̄ = 0.49999999999999
    # ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 10
    D = 10
    I = D^3/10
    EI = E*I

    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
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


    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    kᵤₚ = zeros(2*nᵤ,nₚ)
    kₚₚ = zeros(nₚ,nₚ)
    f = zeros(2*nᵤ)

    opsᵈ[1](elements["Ω"],kᵤᵤ)
    ops[2](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
    ops[3](elements["Ωᵖ"],kₚₚ)
    ops[5](elements["Γᵍ"],kᵤᵤ,f)
    ops[4](elements["Γᵗ"],f)
    k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
    f = [f;zeros(nₚ)]

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]

    push!(nodes,:d₁=>d₁,:d₂=>d₂)

    prescribe!(elements["Ω"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Ω"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
    prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
    prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
    prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)

    h1,l2 = ops[6](elements["Ω"])
    L2 = log10(l2)
    H1 = log10(h1)
    h = ndiv/ndiv_p

    index = 10:40
    XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
        Sheet = xf[1]
        ind = findfirst(n->n==ndiv_p,index)+1
        Sheet["B"*string(ind)] = h
        Sheet["C"*string(ind)] = L2
        Sheet["D"*string(ind)] = H1
    end
end
