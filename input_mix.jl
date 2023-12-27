using Statistics, DelimitedFiles


function import_beam(filename::String)
    elms,~ = ApproxOperator.importmsh(filename1)
    nₚ = length(elms["Ω"][1].x)
    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z

    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]

    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    parameters = (:Linear1D,:□,:CubicSpline)
    n𝒑 = 3

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_NN = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Piecewise{:Seg2},:SegGI2,data)
    f_vN = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegGI2,data)
    f_Nv = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegGI2,data)
    f_MM = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegGI2,data)
    f_vM = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegGI2,data)
    f_Mv = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegGI2,data)
    f_vu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:PoiGI1,data)
    f_Nu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:PoiGI1,data)
    f_Mu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:PoiGI1,data)
    f_vθ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:PoiGI1,data)
    f_Mθ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:PoiGI1,data)
    f_vt = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Poi1},:PoiGI1,data)
    f_vm = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Poi1},:PoiGI1,data)
    f_vb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegRK3,data)
    f_Nb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegRK3,data)
    f_Mb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Seg2},:SegRK3,data)

    elms["Γᵛ∩Ω"] = elms["Γᵛ"]∩elms["Ω"]
    elms["Γᶿ∩Ω"] = elms["Γᶿ"]∩elms["Ω"]
    elms["Γᴹ∩Ω"] = elms["Γᴹ"]∩elms["Ω"]
    elms["Γᵗ∩Ω"] = elms["Γᵗ"]∩elms["Ω"]
    elements["Ω"] = f_Ω(elms["Ω"],sp)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])

    for ap in elements["Ω_Nv"]
        𝓖 = ap.𝓖
        for ξ in 𝓖

    return elements, nodes
end

function import_quad_PP(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
    nₚ = length(elms["Ω"][1].x)
    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z

    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4,data)
    f_Ωᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI16,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant2D,:Quad},:QuadGI4,data)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵍ"] = f_Ωᵍ(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    ) 
    push!(f_Ωᵍ,
    :𝝭=>:𝑠,
    :∂𝝭∂x=>:𝑠,
    :∂𝝭∂y=>:𝑠,
        )
    push!(f_Ωᵖ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠,
    )
    if haskey(elms,"Γᵗ")
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        n₁ = zeros(length(elms["Γᵗ"]))
        n₂ = zeros(length(elms["Γᵗ"]))
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
            :n₁=>(:𝐶,n₁),
            :n₂=>(:𝐶,n₂),
        )
        for ap in elements["Γᵗ"]
            nd₁,nd₂ = ap.𝓒
            x₁ = nd₁.x
            x₂ = nd₂.x
            y₁ = nd₁.y
            y₂ = nd₂.y
            𝐿 = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
            ap.n₁ = (y₂-y₁)/𝐿
            ap.n₂ = (x₁-x₂)/𝐿
        end
    end
    return elements, nodes
end