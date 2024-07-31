
using Statistics , DelimitedFiles
import Gmsh: gmsh

function import_cantilever_Q4P1(filename::String)
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
    # type = PiecewiseParametric{:Constant2D}
    elements["Ωᵖ"] = getPiecewiseElements( entities["Ω"], type, integrationOrder_Ω;)
    elements["Ωᵍᵖ"] = getPiecewiseElements( entities["Ω"], type,  integrationOrder_Ωᵍ;)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    
    return elements, nodes
end
function import_cantilever_Q8P3(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 4
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
    type = PiecewisePolynomial{:Linear2D}
    elements["Ωᵖ"] = getPiecewiseElements( entities["Ω"], type, integrationOrder_Ω )
    elements["Ωᵍᵖ"] =getPiecewiseElements( entities["Ω"], type,  integrationOrder_Ωᵍ)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    
    # gmsh.finalize()
    return elements, nodes
end
function import_cantilever_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 6
    integrationOrder_Γ = 6
    integrationOrder_Ωᵍ =10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    # elements["Γ₁"] = getElements(nodes, entities["Γ₁"],   integrationOrder_Γ)
    # elements["Γ₃"] = getElements(nodes, entities["Γ₃"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Γ₁"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Γ₃"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s =1.8*12/ndiv_p*ones(length(nodes_p))
    # s = 1.3/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
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
    return elements, nodes, nodes_p,Ω,xᵖ,yᵖ,zᵖ, sp,type
end
function import_cantilever_reduce(filename1::String)
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
    # elements["Γ₁"] = getElements(nodes, entities["Γ₁"],   integrationOrder_Γ)
    # elements["Γ₃"] = getElements(nodes, entities["Γ₃"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Γ₁"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Γ₃"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    integrationOrder_1 = 0
    elements["Ωᵛ"] =getElements(nodes, entities["Ω"],   integrationOrder_1)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
    push!(elements["Ωᵛ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    # gmsh.finalize()
    return elements, nodes, nodes_p
end

function import_cantilever_mix_HR(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    integrationOrder_Ω = 10
    integrationOrder_Γ = 10
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
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationOrder_Γ, normal = true)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"], integrationOrder_Γ, normal = true)
    
    
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    
    # type = PiecewisePolynomial{:Constant2D}
    type = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    elements["Γˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭=>:𝑠)
    push!(elements["∂Ωˢ"], :𝝭=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵖ"] = getElements(nodes_p, entities["Γᵍ"], type, integrationOrder_Γ, sp)
    
    nₘ = 6
    # nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Γᵖ"], :𝝭=>:𝑠)
    push!(elements["Γᵖ"], :𝗠=>𝗠)
   

    return elements, nodes, nodes_p, Ω
end

function import_cantilever_mix_bubble(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    integrationOrder_Ω = 8
    integrationOrder_Γ = 8
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
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationOrder_Γ, normal = true)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"], integrationOrder_Γ, normal = true)
    
    
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    
    # # type = PiecewisePolynomial{:Constant2D}
    # type = PiecewisePolynomial{:Linear2D}
    # elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    # elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    # elements["Γˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    # push!(elements["Ωˢ"], :𝝭=>:𝑠)
    # push!(elements["∂Ωˢ"], :𝝭=>:𝑠)

    
    # type = PiecewiseParametric{:Bubble,:Tri3}
      type = PiecewiseParametric{:Bubble,:Quad}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵖ"] = getElements(nodes_p, entities["Γᵍ"], type, integrationOrder_Γ, sp)
    
    nₘ = 6
    # nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Γᵖ"], :𝝭=>:𝑠)
    push!(elements["Γᵖ"], :𝗠=>𝗠)
   

    return elements, nodes, nodes_p, Ω
end

function import_cantilever_mix_internal(filename1::String,filename2::String,filename3::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 5
    integrationOrder_Γ = 5
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
    s = 2.5*s*ones(length(nodes_p))
    # s =1.8*12/ndiv_p*ones(length(nodes_p))
    # s = 1.3/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    gmsh.open(filename3)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
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
function import_cantilever_mix_LM(filename1::String,filename2::String,filename3::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 5
    integrationOrder_Γ = 5
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
    gmsh.open(filename3)
    nodes_λ = get𝑿ᵢ()
    x_λ = nodes_λ.x
    y_λ = nodes_λ.y
    z_λ = nodes_λ.z
    elements["Γ_λ"] = getElements(nodes_λ, entities["Γᵍ"],   integrationOrder_Γ)
    push!(elements["Γ_λ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s =1.8*12/ndiv_p*ones(length(nodes_p))
    # s = 1.3/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)

    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
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
    return elements, nodes, nodes_p,Ω,xᵖ,yᵖ,zᵖ, sp,type
end
function import_cantilever_mix_LM_internal(filename1::String,filename2::String,filename3::String,filename4::String)
    gmsh.initialize()
    gmsh.open(filename1)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    integrationOrder_Ω = 5
    integrationOrder_Γ = 5
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
    gmsh.open(filename3)
    nodes_λ = get𝑿ᵢ()
    x_λ = nodes_λ.x
    y_λ = nodes_λ.y
    z_λ = nodes_λ.z
    elements["Γ_λ"] = getElements(nodes_λ, entities["Γᵍ"],   integrationOrder_Γ)
    push!(elements["Γ_λ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    # s =1.8*12/ndiv_p*ones(length(nodes_p))
    # s = 1.3/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)

    gmsh.open(filename4)
 
    entities = getPhysicalGroups()
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ω, sp)
    # elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"], type,  integrationOrder_Γ, sp)
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

function import_cantilever_T6P3(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)
 
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    # sp = RegularGrid(x,y,z,n = 1,γ = 5)
    integrationOrder_Ω = 5
    integrationOrder_Γ = 5
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
    # gmsh.finalize()
    return elements, nodes, nodes_p
end
function import_cantilever_fem(filename::String)
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
    elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder_Ω)
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder_Γ)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder_Γ)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"],   integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    
    # type = PiecewiseParametric{:Bubble,:Tri3}
      type = PiecewiseParametric{:Bubble,:Quad}
    elements["Ωᵇ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    push!(elements["Ωᵇ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end
prescribeForGauss = quote
   
    
    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*y)
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*y*ν)

end

    

prescribeForPenalty = quote


 
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2)) 
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    # prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍᵖ"],:p₁=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:p₂=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:n₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γᵍᵖ"],:n₂=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

end

prescribeForSquare = quote

    prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->-P*(y-L/2)/6/EI*((6*L-3x)*x + (2+ν)*((y-L/2)^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->P/6/EI*(3*ν*(y-L/2)^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->-P/EI*(L-x)*(y-L/2))
    prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*(y-L/2)^2-D^2/4)))
    prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*(y-L/2)^2 + (4+5*ν)*D^2/4))
    prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->P/EI*(L-x)*(y-L/2)*ν)


    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-(y-L/2)^2)) 
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*(y-L/2)/6/EI*((6*L-3x)*x + (2+ν)*((y-L/2)^2-D^2/4)))
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*(y-L/2)^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    # prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍᵖ"],:p₁=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:p₂=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:n₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γᵍᵖ"],:n₂=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    
end



prescribeForDisplacement = quote
    prescribe!(elements["Γᵗ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γᵗ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2)) 
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    prescribe!(elements["Γ₁"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γ₁"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2)) 
    prescribe!(elements["Γ₃"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    prescribe!(elements["Γ₃"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    # prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γᵍᵖ"],:p₁=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:p₂=>(x,y,z)->-P/EI*(L-x)*y/2)
    # prescribe!(elements["Γᵍᵖ"],:n₁=>(x,y,z)->1.0)
    # prescribe!(elements["Γᵍᵖ"],:n₂=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    prescribe!(elements["Γ₁"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ₁"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ₁"],:n₂₂=>(x,y,z)->1.0)

    prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z)->1.0)

    prescribe!(elements["Γ₃"],:n₁₁=>(x,y,z)->1.0)
    prescribe!(elements["Γ₃"],:n₁₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ₃"],:n₂₂=>(x,y,z)->1.0)




    
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