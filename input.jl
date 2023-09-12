
function import_mf_tri3(filename1::String,filename2::String)
    elms, = ApproxOperator.importmsh(filename1)
    elms_p, = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z

    nodes = Node{(:𝐼,),1}[]
   
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]
    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    # scheme = ApproxOperator.quadraturerule(s)

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2)

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
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        )
    end
    # return elements, nodes, nodes_p
    return elements, nodes
end

function import_fem_tri3(filename1::String,filename2::String)
    elms,nds = ApproxOperator.importmsh(filename1)
    elms_p,nds_p = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    # scheme = ApproxOperator.quadraturerule(s)

    d = zeros(nₚ)
    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI3)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωᵖ,
        :d=>(:𝐼,d),
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
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
        )
    end
    return elements, d
end

function import_quad(filename1::String,filename2::String)
    elms,nds = ApproxOperator.importmsh(filename1)
    elms_p,nds_p = ApproxOperator.importmsh(filename2)
    nₚ = length(elms["Ω"][1].x)

    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    # sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)
    # parameters = (:Linear2D,:□,:CubicSpline)
    # n𝒑 = 21
    # # scheme = ApproxOperator.quadraturerule(s)

    d = zeros(nₚ)

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωᵖ,
        :d=>(:𝐼,d),
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Γᵍ,
        :𝝭=>:𝑠,
    )
    if haskey(elms,"Γᵗ")
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
        )
    end
    return elements, d
end