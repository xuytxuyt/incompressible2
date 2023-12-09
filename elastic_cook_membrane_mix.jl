
using ApproxOperator, LinearAlgebra, Printf ,XLSX
include("input.jl")
# for i in 10:30
ndiv= 20
ndiv_p=19
# elements, nodes, nodes_𝑝,elms = import_rkgsi_mix_quadratic(fid_𝑢,fid_𝑝)
elements,nodes,nodes_p = import_fem_tri3("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(ndiv_p)*".msh")


nᵤ = length(nodes)
nₚ = length(nodes_p)

s = 4.5*10/ndiv_p*ones(nₚ)
push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)



# κ = 400942
# μ = 80.1938
# E = 9*κ*μ/(3*κ+μ)
# ν = (3*κ-2*μ)/2/(3*κ+μ)
E = 70.0
ν = 0.3333
# ν =0.499999


set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Ωᵍ"])
# set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])



ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->6.25)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫qpdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e9*E),
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
    # q = -inv(kₚₚ)*kᵤₚ'*d
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    # push!(nodes,:d₁=>d₁,:d₂=>d₂,:q=>q)




# # fo = open("./vtk/cook_membrance_rkgsi_mix_"*string(ndiv_𝑢)*".vtk","w")
# fo = open("./vtk/cook_membrance_rkgsi_"*string(ndiv_𝑢)*".vtk","w")
# @printf fo "# vtk DataFile Version 2.0\n"
# @printf fo "cook_membrance_rkgsi_mix\n"
# @printf fo "ASCII\n"
# @printf fo "DATASET POLYDATA\n"
# @printf fo "POINTS %i float\n" nₚ
# for p in nodes
#     @printf fo "%f %f %f\n" p.x p.y p.z
# end
# @printf fo "POLYGONS %i %i\n" nₑ 4*nₑ
# for ap in elms["Ω"]
#     𝓒 = ap.vertices
#     @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
# end
# @printf fo "POINT_DATA %i\n" nₚ
# @printf fo "VECTORS U float\n"
# for p in elements["Ωᶜ"]
#     ξ = collect(p.𝓖)[1]
#     N = ξ[:𝝭]
#     u₁ = 0.0
#     u₂ = 0.0
#     for (i,x) in enumerate(p.𝓒)
#         u₁ += N[i]*x.d₁
#         u₂ += N[i]*x.d₂
#     end
#     @printf fo "%f %f %f\n" u₁ u₂ 0.0
# end

# @printf fo "TENSORS STRESS float\n"
# for p in elements["Ωᶜ"]
#     𝓒 = p.𝓒
#     𝓖 = p.𝓖
#     ε₁₁ = 0.0
#     ε₂₂ = 0.0
#     ε₁₂ = 0.0

#     for (i,ξ) in enumerate(𝓖)
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         for (j,xⱼ) in enumerate(𝓒)
#             ε₁₁ += B₁[j]*xⱼ.d₁
#             ε₂₂ += B₂[j]*xⱼ.d₂
#             ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
#         end
#     end
#     σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
#     σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
#     σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
#     @printf fo "%f %f %f\n" σ₁₁ σ₁₂ 0.0
#     @printf fo "%f %f %f\n" σ₁₂ σ₂₂ 0.0
#     @printf fo "%f %f %f\n" 0.0 0.0 0.0
# end
# close(fo)

a = elements["Ω"][end]
ξs = collect(a.𝓖)
𝝭 = ξs[3][:𝝭]
u₂ = 0.0
for (i,x) in enumerate(a.𝓒)
    global u₂ += 𝝭[i]*x.d₂
end
h = nᵤ/nₚ
println(u₂)
println(h)
# index = 10:30
#     XLSX.openxlsx("./xlsx/cook.xlsx", mode="rw") do xf
#         Sheet = xf[5]
#         ind = findfirst(n->n==ndiv_p,index)+1
#         Sheet["B"*string(ind)] = h
#         Sheet["C"*string(ind)] = u₂
       
#     end
# end