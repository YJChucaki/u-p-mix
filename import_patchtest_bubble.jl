
using Gmsh, Statistics

function import_mix_bubble(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    integrationOrder_Ω = 3
    integrationOrder_Γ = 2
    integrationOrder_Ωᵍ = 10

    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s = 2.5*s*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ, normal = true)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ, normal = true)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ, normal = true)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ, normal = true)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    type = PiecewiseParametric{:BubbleTri}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    push!(elements["Ωˢ"], :𝝭=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"], type, integrationOrder_Γ, sp)
    elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"], type, integrationOrder_Γ, sp)
    elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"], type, integrationOrder_Γ, sp)
    elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"], type, integrationOrder_Γ, sp)
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    nₘ = 6
    # nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ⁴ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ⁴ᵖ"], :𝗠=>𝗠)


    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    return elements, nodes, nodes_p
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4/3^0.5*avg𝐴)^0.5
    return s, var𝐴
end

prescribe = quote
    
    prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
    prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))

    prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z)->1.0)

    
    prescribe!(elements["Γ²"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ²"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ²"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ²"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ²"],:n₂₂=>(x,y,z)->1.0)
    

    prescribe!(elements["Γ³"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ³"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ³"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ³"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ³"],:n₂₂=>(x,y,z)->1.0)
    

    prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->1.0)
   
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

end

prescribe_at = quote
    
    prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z)->0.0)
    prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z)->1.0)

    
    prescribe!(elements["Γ²"],:t₁=>(x,y,z)->P)
    prescribe!(elements["Γ²"],:t₂=>(x,y,z)->0.0)

    prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->0.0)
   
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

end

