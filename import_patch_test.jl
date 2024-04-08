
using Tensors, BenchmarkExample
import Gmsh: gmsh

function import_patchtest_up_mix(filename1::String,filename2::String)
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
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"],  integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"],  integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"],  integrationOrder_Γ)
   
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    gmsh.finalize()
    return elements, nodes, nodes_p,xᵖ,yᵖ,zᵖ, sp,type
end

function import_patchtest_up_mix_quad4(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    # sp = RegularGrid(x,y,z,n = 1,γ = 5)
    integrationOrder_Ω = 2
    integrationOrder_Γ = 1
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"],  integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"],  integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"],  integrationOrder_Γ)
   
    
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
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    gmsh.finalize()
    return elements, nodes, nodes_p,xᵖ,yᵖ,zᵖ, sp,type
end
function import_patchtest_fem(filename::String)
    gmsh.initialize()
    gmsh.open(filename)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 2
    integrationOrder_Γ = 1
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"],  integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"],  integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"],  integrationOrder_Γ)
   
    gmsh.finalize()
    return elements, nodes
end
function import_patchtest_mf(filename::String)
    gmsh.initialize()
    gmsh.open(filename)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 2
    integrationOrder_Γ = 1
    integrationOrder_Ωᵍ =10
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 1,γ = 2)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],  type,  integrationOrder_Ω, sp)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],  type, integrationOrder_Ωᵍ, sp)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], type,integrationOrder_Γ, sp)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], type, integrationOrder_Γ, sp)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], type, integrationOrder_Γ, sp)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], type, integrationOrder_Γ, sp)
   
    gmsh.finalize()
    return elements, nodes
end
prescribeForGauss_mix = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))


    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end
prescribeForGauss = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))


    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    
    
end
prescribeForGauss_mix = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))


    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end
prescribeForPenalty = quote
    
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
   

    push!(elements["Γ¹"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ¹"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ²"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ³"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γ⁴"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end
