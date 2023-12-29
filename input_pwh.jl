
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

    𝗠 = zeros(n𝒑)
    ∂𝗠∂x = zeros(n𝒑)
    ∂𝗠∂y = zeros(n𝒑)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    f_Ω = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data)
    f_Ωᵖ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Tri3},:TriGI13,data_p)
    f_Γᵍ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(ReproducingKernel{parameters...,:Seg2},:SegGI5,data)

    elements["Ω"] = f_Ω(elms["Ω"],sp)
    elements["Ωᵖ"] = f_Ωᵖ(elms_p["Ω"],sp_p)
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
function import_fem_tri3_plate_with_hole(filename1::String,filename2::String)
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

    s = 1.5*5/ndiv_p*ones(nₚ)

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
    f_Γᵍ₁ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
    f_Γᵍ₂ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)

    elements["Ω"] = f_Ω(elms["Ω"])
    elements["Ωᵖ"] = f_Ωᵖ(elms["Ω"],sp)
    elements["Γᵍ₁"] = f_Γᵍ₁(elms["Γᵍ₁"])
    elements["Γᵍ₂"] = f_Γᵍ₂(elms["Γᵍ₂"])
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
    push!(f_Γᵍ₁,
        :𝝭=>:𝑠,
    )
    push!(f_Γᵍ₂,
        :𝝭=>:𝑠,
    )
    if haskey(elms,"Γᵗ₁")
        f_Γᵗ₁ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
        elements["Γᵗ₁"] = f_Γᵗ₁(elms["Γᵗ₁"])
        n₁ = zeros(length(elms["Γᵗ₁"]))
        n₂ = zeros(length(elms["Γᵗ₁"]))
        push!(f_Γᵗ₁,
            :𝝭=>:𝑠,
            :n₁=>(:𝐶,n₁),
            :n₂=>(:𝐶,n₂),
        )
        for ap in elements["Γᵗ₁"]
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
    if haskey(elms,"Γᵗ₂")
        f_Γᵗ₂ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
        elements["Γᵗ₂"] = f_Γᵗ₂(elms["Γᵗ₂"])
        n₁ = zeros(length(elms["Γᵗ₂"]))
        n₂ = zeros(length(elms["Γᵗ₂"]))
        push!(f_Γᵗ₂,
            :𝝭=>:𝑠,
            :n₁=>(:𝐶,n₁),
            :n₂=>(:𝐶,n₂),
        )
        for ap in elements["Γᵗ₂"]
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
    if haskey(elms,"Γᵗ₃")
        f_Γᵗ₃ = ApproxOperator.Field{(:𝐼,),1,(:𝑔,:𝐺,:𝐶,:𝑠),4}(Element{:Seg2},:SegGI2,data)
        elements["Γᵗ₃"] = f_Γᵗ₃(elms["Γᵗ₃"])
        n₁ = zeros(length(elms["Γᵗ₃"]))
        n₂ = zeros(length(elms["Γᵗ₃"]))
        push!(f_Γᵗ₃,
            :𝝭=>:𝑠,
            :n₁=>(:𝐶,n₁),
            :n₂=>(:𝐶,n₂),
        )
        for ap in elements["Γᵗ₃"]
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
