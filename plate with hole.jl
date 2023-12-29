
using  ApproxOperator, LinearAlgebra, Printf, XLSX

include("input.jl")

# for i in 1637:1650
    # i=60
    ndiv= 6
    ndiv_p= 6
    # elements,nodes,nodes_p = import_quad("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
    elements,nodes,nodes_p =import_fem_tri3_plate_with_hole("./msh/plate_with_hole_"*string(ndiv)*".msh","./msh/plate_with_hole_"*string(ndiv_p)*".msh")

    nᵤ = length(nodes)
    nₚ = length(nodes_p)

    # s = 1.5*12/ndiv_p*ones(nₚ)

    # push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    # set𝝭!(elements["Ωᵍ"])
    # set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵗ₁"])
    set𝝭!(elements["Γᵗ₂"])
    set𝝭!(elements["Γᵗ₃"])
    set𝝭!(elements["Γᵍ₁"])
    set𝝭!(elements["Γᵍ₂"])

  
    Ē = 3e6
    # ν̄ = 0.49999999999999
    ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
   
  

  T = 1000
  a = 1.0

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))
ApproxOperator.prescribe!(elements["Γᵗ₁"],:t₁=>(x,y,z)->σ₁₁(x,y))
ApproxOperator.prescribe!(elements["Γᵗ₁"],:t₂=>(x,y,z)->σ₁₂(x,y))
ApproxOperator.prescribe!(elements["Γᵗ₂"],:t₁=>(x,y,z)->σ₁₂(x,y))
ApproxOperator.prescribe!(elements["Γᵗ₂"],:t₂=>(x,y,z)->σ₂₂(x,y))
ApproxOperator.prescribe!(elements["Γᵗ₃"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
ApproxOperator.prescribe!(elements["Γᵗ₃"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₂₂=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₁₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₁"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₁"],:g₂=>(x,y,z)->0.0)

ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₂₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₂"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ₂"],:g₂=>(x,y,z)->0.0)
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
    Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄ ),
    ]
    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    kᵤₚ = zeros(2*nᵤ,nₚ)
    kₚₚ = zeros(nₚ,nₚ)
    f = zeros(2*nᵤ)

    ops[3](elements["Ω"],kᵤᵤ)
    ops[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
    ops[5](elements["Ωᵖ"],kₚₚ)
    ops[6](elements["Γᵗ₁"],f)
    ops[6](elements["Γᵗ₂"],f)
    ops[6](elements["Γᵗ₃"],f)
    ops[7](elements["Γᵍ₁"],kᵤᵤ,f)
    ops[7](elements["Γᵍ₂"],kᵤᵤ,f)
    

    k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
    f = [f;zeros(nₚ)]

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    q  = d[2*nᵤ+1:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)

    prescribe!(elements["Ω"],:u=>(x,y,z)->T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) ))
    prescribe!(elements["Ω"],:v=>(x,y,z)->T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) ))
    prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
    

    h1,l2 = ops[10](elements["Ω"],elements["Ωᵖ"])
    # h1,l2 = ops[8](elements["Ω"])
    L2 = log10(l2)
    H1 = log10(h1)
    h = 2nᵤ/nₚ
    println(L2,H1)

#     index = 1637:1650
#     XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
#         Sheet = xf[6]
#         ind = findfirst(n->n==i,index)+1
#         Sheet["B"*string(ind)] = h
#         Sheet["C"*string(ind)] = L2
#         Sheet["D"*string(ind)] = H1
#     end
# end
