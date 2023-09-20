
function import_mf_tri3(filename1::String,filename2::String)
    elms,~= ApproxOperator.importmsh(filename1)
    elms_p,~ = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)
    nᵖ = length(elms_p["Ω"][1].x)

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    xᵖ = elms_p["Ω"][1].x
    yᵖ = elms_p["Ω"][1].y
    zᵖ = elms_p["Ω"][1].z

    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]
    data_p = Dict([:x=>(1,xᵖ),:y=>(1,yᵖ),:z=>(1,zᵖ)])
    nodes_p = [Node{(:𝐼,),1}((i,),data_p) for i in 1:nᵖ]

    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )
    push!(f_Ωᵖ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
    )
    if haskey(elms,"Γᵗ")
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        n₁ = zeros(length(elms["Γᵗ"]))
        n₂ = zeros(length(elms["Γᵗ"]))
        push!(f_Γᵗ,
             :𝝭=>:𝑠,
             :n₁=>(:𝐶,n₁),
             :n₂=>(:𝐶,n₂),     
              :𝗠=>(:𝐶,𝗠),
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

    return elements, nodes, nodes_p
end

function import_fem_tri3(filename1::String,filename2::String)
    elms,~ = ApproxOperator.importmsh(filename1)
    elms_p,~ = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)
    nᵖ = length(elms_p["Ω"][1].x)

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    xᵖ = elms_p["Ω"][1].x
    yᵖ = elms_p["Ω"][1].y
    zᵖ = elms_p["Ω"][1].z
   
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]
    data_p = Dict([:x=>(1,xᵖ),:y=>(1,yᵖ),:z=>(1,zᵖ)])
    nodes_p = [Node{(:𝐼,),1}((i,),data_p) for i in 1:nᵖ]

    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21


    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI3,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωᵖ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
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
    return elements, nodes, nodes_p
end

function import_quad(filename1::String,filename2::String)
    elms,~ = ApproxOperator.importmsh(filename1)
    elms_p,~ = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)
    nᵖ = length(elms_p["Ω"][1].x)
    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    xᵖ = elms_p["Ω"][1].x
    yᵖ = elms_p["Ω"][1].y
    zᵖ = elms_p["Ω"][1].z

    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]
    data_p = Dict([:x=>(1,xᵖ),:y=>(1,yᵖ),:z=>(1,zᵖ)])
    nodes_p = [Node{(:𝐼,),1}((i,),data_p) for i in 1:nᵖ]

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωᵖ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
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
    return elements, nodes, nodes_p
end