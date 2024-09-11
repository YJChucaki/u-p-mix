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
    s = 1.5*s*ones(length(nodes_u))
    push!(nodes_u,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 6

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
    elements["Γᵖ"] = elements["Γ¹ᵗᵖ"]∪elements["Γ²ᵗᵖ"]∪elements["Γ¹ᵍᵖ"]∪elements["Γ²ᵍᵖ"]∪elements["Γ³ᵍᵖ"]
    
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
    elements["Γᵘ"] = elements["Γ¹ᵗᵘ"]∪elements["Γ²ᵗᵘ"]∪elements["Γ¹ᵍᵘ"]∪elements["Γ²ᵍᵘ"]∪elements["Γ³ᵍᵘ"]
   
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
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    # elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γ¹ᵗᵘ"] = getElements(nodes, entities["Γᵗ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵗᵘ"] = getElements(nodes, entities["Γᵗ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ¹ᵍᵘ"] = getElements(nodes, entities["Γᵍ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵍᵘ"] = getElements(nodes, entities["Γᵍ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ³ᵍᵘ"] = getElements(nodes, entities["Γᵍ₃"],  integrationOrder_Γ, normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]
    
    push!(elements["Ωᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["∂Ωᵘ"], :𝝭=>:𝑠)

    push!(elements["Γ¹ᵗᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍᵘ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍᵘ"], :𝝭=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    # elements["∂Ωᵘ"] = getElements(nodes_u, entities["Γ"], type, integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵗᵖ"] = getElements(nodes_p, entities["Γᵗ₁"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵗᵖ"] = getElements(nodes_p, entities["Γᵗ₂"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ¹ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₁"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₂"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵍᵖ"] = getElements(nodes_p, entities["Γᵍ₃"], type, integrationOrder_Γ, sp, normal = true)
    # elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

   
    nₘ = 21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    # push!(elements["∂Ωᵖ"], :𝝭=>:𝑠)
    # push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵗᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍᵖ"], :𝝭=>:𝑠)
    
    push!(elements["Γ¹ᵗᵖ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵗᵖ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵍᵖ"], :𝗠=>𝗠)
   
    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    return elements, nodes, nodes_p , Ω
end

function import_patchtest_Q4P1(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 10
    integrationOrder_Γ = 10
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],  integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γ¹ᵗ"] = getElements(nodes, entities["Γᵗ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵗ"] = getElements(nodes, entities["Γᵗ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ¹ᵍ"] = getElements(nodes, entities["Γᵍ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵍ"] = getElements(nodes, entities["Γᵍ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ³ᵍ"] = getElements(nodes, entities["Γᵍ₃"],  integrationOrder_Γ, normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]
    
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    push!(elements["Γ¹ᵗ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍ"], :𝝭=>:𝑠)
    
    # gmsh.finalize()
    return elements, nodes
end

function import_patchtest_fem(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 10
    integrationOrder_Γ = 10
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],  integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γ¹ᵗ"] = getElements(nodes, entities["Γᵗ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵗ"] = getElements(nodes, entities["Γᵗ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ¹ᵍ"] = getElements(nodes, entities["Γᵍ₁"],  integrationOrder_Γ, normal = true)
    elements["Γ²ᵍ"] = getElements(nodes, entities["Γᵍ₂"],  integrationOrder_Γ, normal = true)
    elements["Γ³ᵍ"] = getElements(nodes, entities["Γᵍ₃"],  integrationOrder_Γ, normal = true)
    # elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]
    
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    push!(elements["Γ¹ᵗ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵗ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵍ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵍ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵍ"], :𝝭=>:𝑠)
    
    # gmsh.finalize()
    return elements, nodes
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
prescribeForFem = quote
    prescribe!(elements["Ω"],:b=>(x,y,z)->b(x,y))
    prescribe!(elements["Γ¹ᵍ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵍ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ³ᵍ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ¹ᵗ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵗ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->T(x,y))
    prescribe!(elements["Ωᵍ"],:p₁=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Ωᵍ"],:p₂=>(x,y,z)->P₂(x,y))
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->P₂(x,y))

    prescribe!(elements["Γ¹ᵗ"],:t=>(x,y,z)->P₂(x,y))
    prescribe!(elements["Γ²ᵍ"],:t=>(x,y,z)->P₁(x,y))
end
prescribe = quote
    
    prescribe!(elements["Ωᵘ"],:b=>(x,y,z)->b(x,y))
    prescribe!(elements["Ωᵖ"],:b=>(x,y,z)->b(x,y))

    prescribe!(elements["Γ¹ᵍᵘ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵍᵘ"],:g=>(x,y,z)->T(x,y))

    prescribe!(elements["Γ³ᵍᵘ"],:g=>(x,y,z)->T(x,y))
    
    prescribe!(elements["Γ¹ᵗᵘ"],:g=>(x,y,z)->T(x,y))
    prescribe!(elements["Γ²ᵗᵘ"],:g=>(x,y,z)->T(x,y))

    prescribe!(elements["Γ¹ᵗᵘ"],:t=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Γ²ᵗᵘ"],:g=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Γ¹ᵗᵘ"],:t=>(x,y,z)->P₂(x,y))
    prescribe!(elements["Γ²ᵗᵘ"],:t₁=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Γ²ᵗᵘ"],:t₂=>(x,y,z)->P₂(x,y))
    
    prescribe!(elements["Γ²ᵍᵘ"],:t=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Γ²ᵍᵘ"],:t=>(x,y,z)->P₂(x,y))
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
    prescribe!(elements["Ωᵍᵖ"],:p₁=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Ωᵍᵖ"],:p₂=>(x,y,z)->P₂(x,y))
    prescribe!(elements["Ωᵍᵖ"],:u=>(x,y,z)->P₁(x,y))
    prescribe!(elements["Ωᵍᵖ"],:v=>(x,y,z)->P₂(x,y))

end

prescribe_Ωᵍᵖ = quote
    prescribe!(elements["Ωᵍᵖ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍᵖ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
end