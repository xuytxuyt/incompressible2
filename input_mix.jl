using Statistics, DelimitedFiles


function import_beam(filename::String)
    elms,~ = ApproxOperator.importmsh(filename)
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

    f_NN = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegGI2,data)
    f_Nv = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegGI2,data)
    f_Nu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Poi1},:PoiGI1,data)
    f_Nb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegRK3,data)
    f_MM = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegGI2,data)
    f_Mv = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegGI2,data)
    f_Mu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Poi1},:PoiGI1,data)
    f_Mb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Seg2},:SegRK3,data)
    f_Mθ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewiseParametric{:Constant1D,:Poi1},:PoiGI1,data)
    f_vN = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data)
    f_vM = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data)
    f_vu = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Poi1},:PoiGI1,data) 
    f_vθ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Poi1},:PoiGI1,data) 
    f_vt = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Poi1},:PoiGI1,data)
    f_vm = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Poi1},:PoiGI1,data)
    f_vb = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegRK3,data)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Poi1},:PoiGI1,data)

    elms["Γᵛ∩Ω"] = elms["Γᵛ"]∩elms["Ω"]
    elms["Γᶿ∩Ω"] = elms["Γᶿ"]∩elms["Ω"]
    elms["Γᴹ∩Ω"] = elms["Γᴹ"]∩elms["Ω"]
    elms["Γᵗ∩Ω"] = elms["Γᵗ"]∩elms["Ω"]

    elements["Ω_NN"] = f_NN(elms["Ω"])
    elements["Ω_Nv"] = f_Nv(elms["Ω"])
    elements["Ω_Nu"] = f_Nu(elms["Γᵛ∩Ω"])
    elements["Ω_Nb"] = f_Nb(elms["Ω"])
    elements["Ω_MM"] = f_MM(elms["Ω"])
    elements["Ω_Mv"] = f_Mv(elms["Ω"])
    # elements["Ω_Mu"] = f_Mu(elms["Γᵛ∩Ω"])
    elements["Ω_Mb"] = f_Mb(elms["Ω"])
    # elements["Ω_Mθ"] = f_Mθ(elms["Γᶿ∩Ω"])
    elements["Ω_vN"] = f_vN(elms["Ω"],sp)
    elements["Ω_vM"] = f_vM(elms["Ω"],sp)
    # elements["Ω_vu"] = f_vu(elms["Γᵛ∩Ω"],sp)
    # elements["Ω_vθ"] = f_vθ(elms["Γᶿ∩Ω"],sp)
    # elements["Ω_vt"] = f_vt(elms["Γᵗ∩Ω"],sp)
    # elements["Ω_vm"] = f_vM(elms["Γᴹ∩Ω"],sp)
    elements["Ω_vb"] = f_vb(elms["Ω"],sp)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])

    push!(f_NN,
        :𝝭=>:𝑠
    )    
    push!(f_Nv,
        :𝝭=>:𝑠
    )
    push!(f_Nu,
        :𝝭=>:𝑠
    )
    push!(f_Nb,
        :𝝭=>:𝑠
    )
    push!(f_MM,
        :𝝭=>:𝑠
    )
    push!(f_Mv,
        :𝝭=>:𝑠
    )
    push!(f_Mu,
        :𝝭=>:𝑠
    )
    push!(f_Mb,
        :𝝭=>:𝑠
    )   
    push!(f_Mθ,
        :𝝭=>:𝑠
    )    
    push!(f_vN,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )        
    push!(f_vM,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )   
    push!(f_vθ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )    
    push!(f_vt,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )    
    push!(f_vm,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )    
    push!(f_vb,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠
    )    
    # for ap in elements["Ω_Nv"]
    #     𝓖 = ap.𝓖
    #     for ξ in 𝓖
    #     end
    # end      
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