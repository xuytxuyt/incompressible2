using Revise, ApproxOperator, LinearAlgebra, Printf, TimerOutputs, XLSX
include("input.jl")

for i in 2:50
    ndiv= i
    ndiv_p= 8
    elements,nodes,nodes_p = import_fem_tri3("./msh/square_"*string(ndiv)*".msh","./msh/square_"*string(ndiv_p)*".msh")
    # elements,nodes,nodes_p = import_fem_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
    nᵤ = length(nodes)
    nₚ = length(nodes_p)

    s = 2.5*10/ndiv_p*ones(nₚ)
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

    E = 3e6
    ν=0.3
    # ν=0.49999999999999
    Ē = E/(1-ν^2)
    ν̄ = ν/(1-ν)

    u(x,y) =  2*x*y+x^2+y^2
    v(x,y) = -2*x*y-x^2-y^2
    ∂u∂x(x,y) = 2*x+2*y
    ∂u∂y(x,y) = 2*x+2*y
    ∂v∂x(x,y) = -2*x-2*y
    ∂v∂y(x,y) = -2*x-2*y
    ∂²u∂x²(x,y) = 2.0
    ∂²u∂x∂y(x,y) = 2.0
    ∂²u∂y²(x,y) = 2.0
    ∂²v∂x²(x,y) = -2.0
    ∂²v∂x∂y(x,y) = -2.0
    ∂²v∂y²(x,y) = -2.0

    ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
    ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
    ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)/(1-2ν)*((1-ν)*∂u∂x(x,y) + ν*∂v∂y(x,y))*n₁+E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₂)
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₁+E/(1+ν)/(1-2ν)*(ν*∂u∂x(x,y) + (1-ν)*∂v∂y(x,y))*n₂)
    ApproxOperator.prescribe!(elements["Ω"],:b₁=>(x,y,z)->-E/(1+ν)/(1-2ν)*((1-ν)*∂²u∂x²(x,y) + ν*∂²v∂x∂y(x,y)) - E/(1+ν)/2*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y)))
    ApproxOperator.prescribe!(elements["Ω"],:b₂=>(x,y,z)->-E/(1+ν)/2*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y)) - E/(1+ν)/(1-2ν)*(ν*∂²u∂x∂y(x,y) + (1-ν)*∂²v∂y²(x,y)))


    ops = [
           Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
           Operator{:∫∫vᵢbᵢdxdy}(),
           Operator{:∫∫p∇vdxdy}(),
           Operator{:∫∫qpdxdy}(:E=>E,:ν=>ν),
           Operator{:∫vᵢtᵢds}(),
           Operator{:∫vᵢgᵢds}(:α=>1e9*E),
           Operator{:Hₑ_PlaneStress}(:E=>Ē,:ν=>ν̄),
    ]
    opsᵛ = [
        Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>E,:ν=>ν )
    ]
    opsᵈ = [
        Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>E,:ν=>ν )
    ]

    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    # kᵤₚ = zeros(2*nᵤ,nₚ)
    # kₚₚ = zeros(nₚ,nₚ)
    kᵤₚ = zeros(2*nᵤ,nᵤ)
    kₚₚ = zeros(nᵤ,nᵤ)
    f = zeros(2*nᵤ)

    opsᵈ[1](elements["Ω"],kᵤᵤ)
    ops[2](elements["Ω"],f)
    # ops[3](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
    ops[3](elements["Ω"],elements["Ω"],kᵤₚ)
    # ops[4](elements["Ωᵖ"],kₚₚ)
    ops[4](elements["Ω"],kₚₚ)  
    ops[5](elements["Γᵗ"],f)
    ops[6](elements["Γᵍ"],kᵤᵤ,f)

    k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
    # f = [f;zeros(nₚ)]
    f = [f;zeros(nᵤ)]

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]

    push!(nodes,:d₁=>d₁,:d₂=>d₂)

    h1,l2 = ops[7](elements["Ω"])
    L2 = log10(l2)
    H1 = log10(h1)
    h = i

    index = 2:50
    XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
        Sheet = xf[7]
        ind = findfirst(n->n==ndiv,index)+1
        Sheet["B"*string(ind)] = h
        Sheet["C"*string(ind)] = L2
        Sheet["D"*string(ind)] = H1
    end
end
