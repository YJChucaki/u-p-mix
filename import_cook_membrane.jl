
using Tensors, BenchmarkExample, Statistics, CairoMakie
import Gmsh: gmsh

function import_cook_membrane_Q4P1(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    integrationOrder_Ω = 1
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 1
    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    type = PiecewisePolynomial{:Constant2D}
    elements["Ωᵖ"] = getMacroElements( entities["Ω"], type, integrationOrder_Ω, 1; )
    elements["Ωᵍᵖ"] = getMacroElements( entities["Ω"], type,  integrationOrder_Ωᵍ, 1;)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    
    gmsh.finalize()
    return elements, nodes
end
function import_cook_membrane_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 2
    integrationOrder_Γ = 2
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
   
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s = 1.5/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # gmsh.finalize()
    return elements, nodes, nodes_p,Ω
end
function import_cook_membrane_mix_tri6(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 4
    integrationOrder_Γ = 4
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)

   
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s = 1.5/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    # gmsh.finalize()
    return elements, nodes, nodes_p,xᵖ,yᵖ,zᵖ, sp,type
end
function import_cook_membrane_mix_quad4(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    # sp = RegularGrid(x,y,z,n = 1,γ = 5)
    integrationOrder_Ω = 3
    integrationOrder_Γ = 2
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)

    
    gmsh.open(filename2)
    
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    nᵖ = length(nodes_p)
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z

    # s, var𝐴 = cal_area_support(nodes_p)
    # s = 1.5*s*ones(nᵖ)
    
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 1,γ = 2)
    gmsh.open(filename1)
    integrationOrder= 3
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    gmsh.finalize()
    return elements, nodes, nodes_p,xᵖ,yᵖ,zᵖ, sp,type
end
function import_cook_membrane_T6P3(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    # sp = RegularGrid(x,y,z,n = 1,γ = 5)
    integrationOrder_Ω = 4
    integrationOrder_Γ = 4
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    nᵖ = length(nodes_p)
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    # s, var𝐴 = cal_area_support(nodes_p)
    # s = 1.5*s*ones(nᵖ)
    # gmsh.open(filename1)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"],  integrationOrder_Ω)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],  integrationOrder_Γ)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"],  integrationOrder_Ωᵍ)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    gmsh.finalize()
    return elements, nodes, nodes_p
end
function import_cook_membrane_fem(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 4
    integrationOrder_Γ = 4
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    
  
   
    
    gmsh.finalize()
    return elements, nodes
end


    

prescribeForPenalty = quote
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->6.25)
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍᵖ"],:n₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍᵖ"],:n₂=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    
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