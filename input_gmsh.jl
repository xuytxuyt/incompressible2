
function import_curved_beam(filename::String)
    elms,~= ApproxOperator.importmsh(filename)
    nₚ = length(elms["Ω"][1].x)
    x = elms["Ω"][1].x
    y = elms["Ω"][1].y
    z = elms["Ω"][1].z
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    nodes = [Node{(:𝐼,),1}((i,),data) for i in 1:nₚ]

    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=1)
    n𝒑 = 3
    𝗠 = zeros(n𝒑)
    
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    as = elms["Ω"]
    ## Ωᵥ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Seg2},:SegRK3,data)
    type = getfield(f,:type)
    elements["Ωᵥ"] = type[]
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    ni = 3
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝑔, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            x_[f.𝐺], y_[f.𝐺], z_[f.𝐺] = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Ωᵥ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    # Ωₙₙ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Seg2},:SegGI2,data)
    type = getfield(f,:type)
    elements["Ωₙₙ"] = type[]
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    ni = 2
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝑔, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),

    )
    𝑛𝑝 = ApproxOperator.get𝑛𝑝(type((0,0,Node{(:𝐼,),1}[]),(0,0,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[])))
    for (C,a) in enumerate(as)
        for i in 1:𝑛𝑝
            f.𝐼 += 1
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            x_[f.𝐺], y_[f.𝐺], z_[f.𝐺] = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Ωₙₙ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Ωₘₘ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Seg2},:SegGI2,data)
    type = getfield(f,:type)
    elements["Ωₘₘ"] = type[]
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    ni = 2
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝑔, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    𝑛𝑝 = ApproxOperator.get𝑛𝑝(type((0,0,Node{(:𝐼,),1}[]),(0,0,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[])))
    for (C,a) in enumerate(as)
        for i in 1:𝑛𝑝
            f.𝐼 += 1
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            x_[f.𝐺], y_[f.𝐺], z_[f.𝐺] = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Ωₘₘ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Ωₙᵥ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Seg2},:SegRK3,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Ωₙᵥ"] = type[]
    ni = 3
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝑔, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    𝑛𝑝 = ApproxOperator.get𝑛𝑝(type((0,0,Node{(:𝐼,),1}[]),(0,0,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[])))
    for (C,a) in enumerate(as)
        for i in 1:𝑛𝑝
            f.𝐼 += 1
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            x_[f.𝐺], y_[f.𝐺], z_[f.𝐺] = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Ωₙᵥ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        )        
    end

    # Ωₘᵥ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Seg2},:SegRK3,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Ωₘᵥ"] = type[]
    ni = 3
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝑔, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    𝑛𝑝 = ApproxOperator.get𝑛𝑝(type((0,0,Node{(:𝐼,),1}[]),(0,0,Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[])))
    for (C,a) in enumerate(as)
        for i in 1:𝑛𝑝
            f.𝐼 += 1
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            x_[f.𝐺], y_[f.𝐺], z_[f.𝐺] = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Ωₘᵥ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    ## Γ
    elms["Γ"] = [ApproxOperator.Poi1((i,),x,y,z) for i in 1:nₚ]
    elms["Γ∩Ω"] = elms["Γ"]∩elms["Ω"]

    abs = elms["Γ∩Ω"]
    # Γᵥ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵥ"] = type[]
    ni = 1
    ne = length(abs)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, ξ_),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,(a,b)) in enumerate(abs)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            ((x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]),ξ_[f.𝐺]) = b(a,ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵥ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    # Γₙ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γₙ"] = type[]
    ni = 1
    ne = length(abs)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    n₁ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, ξ_),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
        :n₁=> (:𝐺, n₁)
    )
    for (C,(a,b)) in enumerate(abs)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            ((x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]),ξ_[f.𝐺]) = b(a,ξ)
            if ξ_[f.𝐺] ≈ -1.0
                n₁[f.𝐺] = -1.0
            else ξ_[f.𝐺] ≈ 1.0
                n₁[f.𝐺] = 1.0
            end
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γₙ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Γₘ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γₘ"] = type[]
    ni = 1
    ne = length(abs)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    n₁ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, ξ_),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
        :n₁=> (:𝐺, n₁)
    )
    for (C,(a,b)) in enumerate(abs)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            ((x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]),ξ_[f.𝐺]) = b(a,ξ)
            if ξ_[f.𝐺] ≈ -1.0
                n₁[f.𝐺] = -1.0
            else ξ_[f.𝐺] ≈ 1.0
                n₁[f.𝐺] = 1.0
            end
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γₘ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Γᵍ
    elms["Γᵛ∩Ω"] = elms["Γᵛ"]∩elms["Ω"]
    as = elms["Γᵛ"]
    abs = elms["Γᵛ∩Ω"]
    # Γᵛᵥ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵛᵥ"] = type[]
    ni = 1
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            (x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]) = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵛᵥ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    # Γᵛₙ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵛₙ"] = type[]
    ni = 1
    ne = length(abs)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    n₁ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, ξ_),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
        :n₁=> (:𝐺, n₁)
    )
    for (C,(a,b)) in enumerate(abs)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            ((x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]),ξ_[f.𝐺]) = b(a,ξ)
            if ξ_[f.𝐺] ≈ -1.0
                n₁[f.𝐺] = -1.0
            else ξ_[f.𝐺] ≈ 1.0
                n₁[f.𝐺] = 1.0
            end
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵛₙ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Γᵛₘ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(PiecewisePolynomial{:Linear1D,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵛₘ"] = type[]
    ni = 1
    ne = length(abs)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    n₁ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, ξ_),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
        :n₁=> (:𝐺, n₁)
    )
    for (C,(a,b)) in enumerate(abs)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            ((x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]),ξ_[f.𝐺]) = b(a,ξ)
            if ξ_[f.𝐺] ≈ -1.0
                n₁[f.𝐺] = -1.0
            else ξ_[f.𝐺] ≈ 1.0
                n₁[f.𝐺] = 1.0
            end
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵛₘ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
        ) 
    end

    # Γᶿ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᶿ"] = type[]
    ni = 1
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            (x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]) = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᶿ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    # Γᵗ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵗ"] = type[]
    ni = 1
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            (x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]) = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵗ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    # Γᵐ
    f = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{:Linear1D,:□,:CubicSpline,:Poi1},:PoiGI1,data)
    data𝓖 = getfield(f,:data𝓖)
    weights = data𝓖[:w][2]
    positions = data𝓖[:ξ][2]
    scheme = zip(weights,positions)
    type = getfield(f,:type)
    elements["Γᵐ"] = type[]
    ni = 1
    ne = length(as)
    𝑤_ = zeros(ni*ne)
    x_ = zeros(ni*ne)
    y_ = zeros(ni*ne)
    z_ = zeros(ni*ne)
    ξ_ = zeros(ni*ne)
    push!(f,
        :w => (:𝑔, weights),
        :ξ => (:𝐺, positions),
        :𝑤 => (:𝐺, 𝑤_),
        :x => (:𝐺, x_),
        :y => (:𝐺, y_),
        :z => (:𝐺, z_),
    )
    for (C,a) in enumerate(as)
        indices = Set{Int}()
        for (w,ξ) in scheme
            xᵢ,yᵢ,zᵢ = a(ξ)
            union!(indices,sp(xᵢ,yᵢ,zᵢ))
        end
        ni = length(indices)
        for i in indices
            f.𝐼 = i
            ApproxOperator.add𝓒!(f)
        end
        for (g,(w,ξ)) in enumerate(scheme)
            f.𝑔 = g
            f.𝐺 += 1
            f.𝐶 = C
            ApproxOperator.add𝓖!(f)
            f.𝑠 += ni
            𝑤_[f.𝐺] = ApproxOperator.get𝐽(a,ξ)*w
            (x_[f.𝐺], y_[f.𝐺], z_[f.𝐺]) = a(ξ)
        end
        𝓒 = ApproxOperator.get𝓒(f)
        𝓖 = ApproxOperator.get𝓖(f)
        push!(elements["Γᵐ"],type(𝓒,𝓖))
        push!(f,
            :𝝭=>:𝑠,
            :𝗠=>(:𝐶,𝗠),
        ) 
    end

    return elements,nodes, elms
end