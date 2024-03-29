
using Tensors, BenchmarkExample
import Gmsh: gmsh

function import_cantilever_mix_tri3(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
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
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    gmsh.finalize()
    return elements, nodes, nodes_p,xᵖ,yᵖ,zᵖ, sp,type
end

function import_cantilever_mix_quad4(filename1::String,filename2::String)
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
prescribeForGauss = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)

    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    push!(elements["Ωᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end

prescribeForPenalty = quote
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Γᵍᵖ"],:p₁=>(x,y,z)->-P/EI*(L-x)*y/2)
    prescribe!(elements["Γᵍᵖ"],:p₂=>(x,y,z)->-P/EI*(L-x)*y/2)
    prescribe!(elements["Γᵍᵖ"],:n₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍᵖ"],:n₂=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    push!(elements["Γᵗ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵍ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    
end
