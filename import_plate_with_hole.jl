using Gmsh, Statistics
using CairoMakie
# using GLMakie
function import_patchtest_mix(filename1::String, filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_u = get𝑿ᵢ()
    xᵘ = nodes_u.x
    yᵘ = nodes_u.y
    zᵘ = nodes_u.z
    Ω = getElements(nodes_u, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes_u))
    push!(nodes_u,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 5
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 5

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵖ"] = getElements(nodes, entities["Ω"],  integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γ¹ᵗᵖ"] = getElements(nodes, entities["Γᵗ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵗᵖ"] = getElements(nodes, entities["Γᵗ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ¹ᵍᵖ"] = getElements(nodes, entities["Γᵍ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵍᵖ"] = getElements(nodes, entities["Γᵍ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ³ᵍᵖ"] = getElements(nodes, entities["Γᵍ₃"],  integrationOrder_Γ, normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]
    
    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["∂Ωᵖ"], :𝝭=>:𝑠)

    push!(elements["Γ¹ᵗᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍᵖ"], :𝝭=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵘ,yᵘ,zᵘ,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes_u, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵘ"] = getElements(nodes_u, entities["Γ"], type, integrationOrder_Γ, sp)
    elements["Ωᵍᵘ"] = getElements(nodes_u, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵗᵘ"] = getElements(nodes_u, entities["Γᵗ₁"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵗᵘ"] = getElements(nodes_u, entities["Γᵗ₂"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ¹ᵍᵘ"] = getElements(nodes_u, entities["Γᵍ₁"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵍᵘ"] = getElements(nodes_u, entities["Γᵍ₂"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵍᵘ"] = getElements(nodes_u, entities["Γᵍ₃"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

   
    nₘ = 21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["∂Ωᵘ"], :𝝭=>:𝑠)
    push!(elements["∂Ωᵘ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵗᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍᵘ"], :𝝭=>:𝑠)
    
    push!(elements["Γ¹ᵗᵘ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵗᵘ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵍᵘ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵍᵘ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵍᵘ"], :𝗠=>𝗠)
   
    push!(elements["Ωᵘ"], :𝝭=>:𝑠)
    push!(elements["Ωᵘ"],  :𝗠=>𝗠)
    push!(elements["Ωᵍᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    gmsh.finalize()
    return elements, nodes, nodes_u
end

function import_patchtest_mix_old(filename1::String, filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes_p))
    # s = 1.5/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 10
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 10
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes, entities["Ω"],  integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ, normal = true)
    elements["Γ¹ᵘ"] = getElements(nodes, entities["Γ¹"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵘ"] = getElements(nodes, entities["Γ²"],  integrationOrder_Γ, normal = true)
    elements["Γ³ᵘ"] = getElements(nodes, entities["Γ³"],  integrationOrder_Γ, normal = true)
    elements["Γ⁴ᵘ"] = getElements(nodes, entities["Γ⁴"],  integrationOrder_Γ, normal = true)
    elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

    
    push!(elements["Ωᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ¹ᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ²ᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ³ᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ⁴ᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)


    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp, normal = true)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp, normal = true)
    elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

   
    nₘ = 21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Γ¹ᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ¹ᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ²ᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ²ᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ³ᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ³ᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ⁴ᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ⁴ᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

   
    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # gmsh.finalize()
    return elements, nodes, nodes_p , Ω
end


function resetx!(a::ApproxOperator.AbstractElement)
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    x₁ = 𝓒[1].x
    x₂ = 𝓒[2].x
    x₃ = 𝓒[3].x
    y₁ = 𝓒[1].y
    y₂ = 𝓒[2].y
    y₃ = 𝓒[3].y
    z₁ = 𝓒[1].z
    z₂ = 𝓒[2].z
    z₃ = 𝓒[3].z
    for ξ in 𝓖
        N₁ = ξ.ξ
        N₂ = ξ.η
        N₃ = 1.0-ξ.ξ-ξ.η
        ξ.x = N₁*x₁ + N₂*x₂ + N₃*x₃
        ξ.y = N₁*y₁ + N₂*y₂ + N₃*y₃
        ξ.z = N₁*z₁ + N₂*z₂ + N₃*z₃
    end
end

function reset𝐴!(a::ApproxOperator.AbstractElement)
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    x₁ = 𝓒[1].x
    x₂ = 𝓒[2].x
    x₃ = 𝓒[3].x
    y₁ = 𝓒[1].y
    y₂ = 𝓒[2].y
    y₃ = 𝓒[3].y
    a.𝐴 = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
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
    
    prescribe!(elements["Ωᵘ"],:b=>(x,y,z)->b(x,y))
    prescribe!(elements["Ωᵖ"],:b=>(x,y,z)->b(x,y))

    prescribe!(elements["Γ¹ᵍᵘ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵍᵘ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ³ᵍᵘ"],:g=>(x,y,z)->T(x,y))
    
    prescribe!(elements["Γ¹ᵗᵘ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵗᵘ"],:g=>(x,y,z)->T(x,y))
   
    # prescribe!(elements["Γ¹ᵖ"],:n₁₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γ¹ᵖ"],:n₁₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ¹ᵖ"],:n₂₂=>(x,y,z)->1.0)

    # prescribe!(elements["Γ²ᵖ"],:n₁₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γ²ᵖ"],:n₁₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ²ᵖ"],:n₂₂=>(x,y,z)->1.0)

    # prescribe!(elements["Γ³ᵖ"],:n₁₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γ³ᵖ"],:n₁₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ³ᵖ"],:n₂₂=>(x,y,z)->1.0)

    # prescribe!(elements["Γ⁴ᵖ"],:n₁₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γ⁴ᵖ"],:n₁₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ⁴ᵖ"],:n₂₂=>(x,y,z)->1.0)

    prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->T(x,y))
    prescribe!(elements["Ωᵍᵖ"],:u=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Ωᵍᵖ"],:v=>(x,y,z)->P₂(x,y))

end
prescribeForFem = quote
    
    prescribe!(elements["Ω"],:s=>(x,y,z)->s(x,y))

    prescribe!(elements["Γ¹"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ³"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ⁴"],:g=>(x,y,z)->T(x,y))
   

    
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->T(x,y))

end

prescribe_Ωᵍᵖ = quote
    prescribe!(elements["Ωᵍᵖ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍᵖ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
end