using Statistics, DelimitedFiles

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
    sp_p = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    s, var𝐴 = cal_area_support(elms["Ω"])
    sₚ, var𝐴 = cal_area_support(elms_p["Ω"])
    sₚ= 1.5*sₚ*ones(nᵖ)
    s= 1.5*s*ones(nₚ)

    push!(nodes_p,:s₁=>sₚ,:s₂=>sₚ,:s₃=>sₚ)
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI5,data)

    elements["Ω"] = f_Ω(elms["Ω"],sp)
    elements["Ωᵖ"] = f_Ωᵖ(elms["Ω"],sp_p)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"],sp)
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
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI5,data)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"],sp)
        n₁ = zeros(length(elms["Γᵗ"]))
        n₂ = zeros(length(elms["Γᵗ"]))
        push!(f_Γᵗ,
             :𝝭=>:𝑠,
             :n₁=>(:𝐶,n₁),
             :n₂=>(:𝐶,n₂),     
             :𝗠=>(:𝐶,𝗠),
        )
        for (ap,a) in zip(elements["Γᵗ"],elms["Γᵗ"])
            x₁ = a.x[a.i[1]]
            x₂ = a.x[a.i[2]]
            y₁ = a.y[a.i[1]]
            y₂ = a.y[a.i[2]]
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

    s, var𝐴 = cal_area_support(elms_p["Ω"])
    s = 1.5*s*ones(nᵖ)

    f = open("./xlsx/var.txt", "a")
    writedlm(f, [nᵖ var𝐴])
    
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    sp = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21


    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI13,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms["Ω"],sp)
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

function import_fem_tri3_GI1(filename1::String,filename2::String)
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

    sp = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21


    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI13,data)
    f_Ωᵛ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI1,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵛ"] = f_Ωᵛ(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"],sp)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωᵛ,
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
function import_fem_tri3_direct(filename1::String,filename2::String)
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

    sp = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21


    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Tri3},:TriGI3,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    # f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Poi1},:PoiGI1,data)
    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"],sp)
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
    # push!(f_Γᵍ,
    #     :𝝭=>:𝑠,
    # )

    elements["Γᵍ"] = Element{:Poi1}[]
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    nₑ = length(elms["Γᵍ"])

    for (C,a) in enumerate(elms["Γᵍ"])
        element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
        push!(𝓒,nodes[a.i[1]])
        push!(elements["Γᵍ"],element)
        c += 1
        if C == nₑ
            element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
            push!(𝓒,nodes[a.i[2]])    
            push!(elements["Γᵍ"],element)
        end
    end
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

    sp = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4,data)
    f_Ωᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI16,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Quad},:QuadGI4,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵍ"] = f_Ωᵍ(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms["Ω"],sp)
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

function import_quad_GI1(filename1::String,filename2::String)
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

    sp = ApproxOperator.RegularGrid(xᵖ,yᵖ,zᵖ,n=1,γ=2)
    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI4,data)
    f_Ωᵛ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad},:QuadGI1,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Quad},:QuadGI4,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵛ"] = f_Ωᵛ(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"],sp)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
        push!(f_Ωᵛ,
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

function import_quad8_GI1(filename1::String,filename2::String)
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


    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad8},:QuadGI9,data)
    f_Ωᵛ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad8},:QuadGI1,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Quad8},:QuadGI9,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg3},:SegGI3,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵛ"] = f_Ωᵛ(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"])
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
        push!(f_Ωᵛ,
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
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg3},:SegGI3,data)
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

function import_fem_bar(filename1::String,filename2::String,filename3::String)
    elms,~ = ApproxOperator.importmsh(filename1)
    elms_n,~ = ApproxOperator.importmsh(filename2)
    elms_v,~ = ApproxOperator.importmsh(filename3)
    nₚ = length(elms["Ω"][1].x)
    nₙ = length(elms_n["Ω"][1].x)
    nᵥ = length(elms_v["Ω"][1].x)
    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    xₙ = elms_n["Ω"][1].x
    yₙ = elms_n["Ω"][1].y
    zₙ = elms_n["Ω"][1].z
    xᵥ = elms_v["Ω"][1].x
    yᵥ = elms_v["Ω"][1].y
    zᵥ = elms_v["Ω"][1].z

    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]
    data_n = Dict([:x=>(1,xₙ),:y=>(1,yₙ),:z=>(1,zₙ)])
    nodes_n = [Node{(:𝐼,),1}((i,),data_n) for i in 1:nₙ]
    data_v = Dict([:x=>(1,xᵥ),:y=>(1,yᵥ),:z=>(1,zᵥ)])
    nodes_v = [Node{(:𝐼,),1}((i,),data_v) for i in 1:nᵥ]

    sp_n = ApproxOperator.RegularGrid(xₙ,yₙ,zₙ,n=1,γ=2)
    sp_v = ApproxOperator.RegularGrid(xᵥ,yᵥ,zᵥ,n=1,γ=2)
    parameters = (:Linear1D,:□,:CubicSpline)
    n𝒑 = 21

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
    f_Ωⁿ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data_n)
    f_Ωᵛ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI2,data_v)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Poi1},:PoiGI1,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωⁿ"] = f_Ωⁿ(elms["Ω"],sp_n)
    elements["Ωᵛ"] = f_Ωᵛ(elms["Ω"],sp_v)
    elements["Γᵍ"] = f_Γᵍ(elms["Γᵍ"])
    push!(f_Ω,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
    )
    push!(f_Ωⁿ,
        :𝝭=>:𝑠,
        :∂𝝭∂x=>:𝑠,
        :∂𝝭∂y=>:𝑠,
        :𝗠=>(:𝐶,𝗠),
        :∂𝗠∂x=>(:𝐶,∂𝗠∂x),
        :∂𝗠∂y=>(:𝐶,∂𝗠∂y)
    )
    push!(f_Ωᵛ,
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
        f_Γᵗ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Poi1},:PoiGI1,data)
        elements["Γᵗ"] = f_Γᵗ(elms["Γᵗ"])
        push!(f_Γᵗ,
            :𝝭=>:𝑠,
        )
    end
    return elements, nodes, nodes_n, nodes_v
end

function cal_area_support(elms::Vector{ApproxOperator.Tri3})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.x[elm.i[1]]
        y₁ = elm.y[elm.i[1]]
        x₂ = elm.x[elm.i[2]]
        y₂ = elm.y[elm.i[2]]
        x₃ = elm.x[elm.i[3]]
        y₃ = elm.y[elm.i[3]]
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = 4/3^0.5*avg𝐴
    return s, var𝐴
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