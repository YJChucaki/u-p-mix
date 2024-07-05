using Gmsh, Statistics
using CairoMakie
# using GLMakie

function import_patchtest_Q4P1(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    integrationOrder_Ω = 3
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 2
    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    type = PiecewisePolynomial{:Constant2D}
    # type = PiecewiseParametric{:Constant2D}
    elements["Ωᵖ"] = getPiecewiseElements( entities["Ω"], type, integrationOrder_Ω;)
    elements["Ωᵍᵖ"] = getPiecewiseElements( entities["Ω"], type,  integrationOrder_Ωᵍ;)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
 
    return elements, nodes
end
function import_patchtest_Q4R1(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    integrationOrder_Ω = 3
    integrationOrder_Ωᵖ = 1
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 2
    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    elements["Ωᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    gmsh.finalize()
    return elements, nodes
end
function import_patchtest_Q8P3(filename::String)
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
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    type = PiecewisePolynomial{:Linear2D}
    elements["Ωᵖ"] = getPiecewiseElements( entities["Ω"], type, integrationOrder_Ω )
    elements["Ωᵍᵖ"] =getPiecewiseElements( entities["Ω"], type,  integrationOrder_Ωᵍ)
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    gmsh.finalize()
    return elements, nodes
end
function import_patchtest_mix(filename1::String, filename2::String)
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

    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 6
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)


    nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    # gmsh.finalize()
    return elements, nodes, nodes_p , Ω
end
function import_patchtest_mix_LM(filename1::String, filename2::String,filename3::String)
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
    s = 1.5*s*ones(length(nodes_p))
    # s = 1.5/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 2
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    
    gmsh.open(filename3)
    nodes_λ = get𝑿ᵢ()
    x_λ = nodes_λ.x
    y_λ = nodes_λ.y
    z_λ = nodes_λ.z
    elements["Γ_λ1"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ_λ2"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ_λ3"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ_λ4"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)

    elements["Γ_λ"] = elements["Γ_λ1"]∪elements["Γ_λ2"]∪elements["Γ_λ3"]∪elements["Γ_λ4"]

    push!(elements["Γ_λ1"], :𝝭=>:𝑠)
    push!(elements["Γ_λ2"], :𝝭=>:𝑠)
    push!(elements["Γ_λ3"], :𝝭=>:𝑠)
    push!(elements["Γ_λ4"], :𝝭=>:𝑠)

    nₘ = 6
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    gmsh.finalize()
    return elements, nodes, nodes_p , Ω
end
function import_patchtest_T6P3(filename1::String, filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 4
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = getfield(nodes_p[1],:data)[:x][2]
    yᵖ = getfield(nodes_p[1],:data)[:y][2]
    zᵖ = getfield(nodes_p[1],:data)[:z][2]
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"],   integrationOrder_Ωᵍ)
    nₘ = 6
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    gmsh.finalize()
    return elements, nodes, nodes_p ,xᵖ,yᵖ,zᵖ
end
function import_patchtest_mix_tri6(filename1::String, filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()
    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = getfield(nodes_p[1],:data)[:x][2]
    yᵖ = getfield(nodes_p[1],:data)[:y][2]
    zᵖ = getfield(nodes_p[1],:data)[:z][2]
    # xᵖ = nodes_p.x
    # yᵖ = nodes_p.y
    # zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 4.0*s*ones(length(nodes_p))
    # s = 1.5/10*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 4
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)


    nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)
    # gmsh.finalize()
    return elements, nodes, nodes_p ,Ω
end
function import_patchtest_quad(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ωᵛ"] = getElements(nodes, entities["Ω"], 0)
    elements["Ωᵈ"] = getElements(nodes, entities["Ω"], 3)

    push!(elements["Ωᵛ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵈ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    gmsh.finalize()

    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    xg = getfield(elements["Ωᵛ"][1].𝓖[1],:data)[:x][2]
    yg = getfield(elements["Ωᵛ"][1].𝓖[1],:data)[:y][2]
    zg = getfield(elements["Ωᵛ"][1].𝓖[1],:data)[:z][2]

    lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
    f = Figure(backgroundcolor = :transparent)
    ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)
    lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

    for elm in elements["Ωᵛ"]
        id = [node.𝐼 for node in elm.𝓒]
        lines!(x[id[[1,2,3,4,1]]],y[id[[1,2,3,4,1]]], linewidth = lwm, color = :black)
    end
    scatter!(x,y,marker = :circle, markersize = mso, color = :black)
    scatter!(xg,yg,marker = :cross, markersize = 5, color = :blue)

    return elements, nodes, f
end


function import_patchtest_stripe(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 2
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    gmsh.finalize()

    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    xg = getfield(elements["Ω"][1].𝓖[1],:data)[:x][2]
    yg = getfield(elements["Ω"][1].𝓖[1],:data)[:y][2]
    zg = getfield(elements["Ω"][1].𝓖[1],:data)[:z][2]

    lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
    f = Figure(backgroundcolor = :transparent)
    ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)
    lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

    for elm in elements["Ω"]
        id = [node.𝐼 for node in elm.𝓒]
        lines!(x[id[[1,2,3,1]]],y[id[[1,2,3,1]]], linewidth = lwm, color = :black)
    end
    scatter!(x,y,marker = :circle, markersize = mso, color = :black)
    scatter!(xg,yg,marker = :cross, markersize = 5, color = :blue)

    return elements, nodes, f
end

function import_patchtest_unionJack(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    integrationOrder = 2
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    stripe2unionJack!(elements["Ω"])

    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    gmsh.finalize()

    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    xg = getfield(elements["Ω"][1].𝓖[1],:data)[:x][2]
    yg = getfield(elements["Ω"][1].𝓖[1],:data)[:y][2]
    zg = getfield(elements["Ω"][1].𝓖[1],:data)[:z][2]

    lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
    f = Figure(backgroundcolor = :transparent)
    ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)
    lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

    for elm in elements["Ω"]
        id = [node.𝐼 for node in elm.𝓒]
        lines!(x[id[[1,2,3,1]]],y[id[[1,2,3,1]]], linewidth = lwm, color = :black)
    end
    scatter!(x,y,marker = :circle, markersize = mso, color = :black)
    scatter!(xg,yg,marker = :cross, markersize = 5, color = :blue)

    return elements, nodes, f
end

function import_patchtest_cross(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 2
 
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    stripe2cross!(elements["Ω"], nodes)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)
    gmsh.finalize()

    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    xg = getfield(elements["Ω"][1].𝓖[1],:data)[:x][2]
    yg = getfield(elements["Ω"][1].𝓖[1],:data)[:y][2]
    zg = getfield(elements["Ω"][1].𝓖[1],:data)[:z][2]

    lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
    f = Figure(backgroundcolor = :transparent)
    ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)
    lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

    for elm in elements["Ω"]
        id = [node.𝐼 for node in elm.𝓒]
        lines!(x[id[[1,2,3,1]]],y[id[[1,2,3,1]]], linewidth = lwm, color = :black)
    end
    scatter!(x,y,marker = :circle, markersize = mso, color = :black)
    scatter!(xg,yg,marker = :cross, markersize = 5, color = :blue)

    return elements, nodes, f
end

function import_patchtest_fem(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    integrationOrder_Ω = 6
    integrationOrder_Ωᵍ = 10
    integrationOrder_Γ = 6
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    # gmsh.finalize()

    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    xg = getfield(elements["Ω"][1].𝓖[1],:data)[:x][2]
    yg = getfield(elements["Ω"][1].𝓖[1],:data)[:y][2]
    zg = getfield(elements["Ω"][1].𝓖[1],:data)[:z][2]

    # lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
    # f = Figure(backgroundcolor = :transparent)
    # ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
    # hidespines!(ax)
    # hidedecorations!(ax)
    # lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

    # for elm in elements["Ω"]
    #     id = [node.𝐼 for node in elm.𝓒]
    #     lines!(x[id[[1,2,3,1]]],y[id[[1,2,3,1]]], linewidth = lwm, color = :black)
    # end
    # scatter!(x,y,marker = :circle, markersize = mso, color = :black)
    # scatter!(xg,yg,marker = :cross, markersize = 5, color = :blue)

    return elements, nodes
end

function stripe2unionJack!(elms::Vector{ApproxOperator.AbstractElement})
    nₑ = length(elms)
    nᵢ = Int((nₑ/2)^0.5/2)
    for i in 1:nᵢ
        for j in 1:nᵢ
            ind = 8*nᵢ*(j-1)+4*i-3
            x₁ = elms[ind].𝓒[1]
            x₂ = elms[ind].𝓒[2]
            x₃ = elms[ind+1].𝓒[3]
            x₄ = elms[ind+1].𝓒[1]
            elms[ind].𝓒[1:3] .= [x₁, x₂, x₃]
            resetx!(elms[ind])
            elms[ind+1].𝓒[1:3] .= [x₄, x₁, x₃]
            resetx!(elms[ind+1])

            ind = 8*nᵢ*(j-1)+4*nᵢ+4*i-1
            x₁ = elms[ind].𝓒[1]
            x₂ = elms[ind].𝓒[2]
            x₃ = elms[ind+1].𝓒[3]
            x₄ = elms[ind+1].𝓒[1]
            elms[ind].𝓒[1:3] .= [x₁, x₂, x₃]
            resetx!(elms[ind])
            elms[ind+1].𝓒[1:3] .= [x₄, x₁, x₃]
            resetx!(elms[ind+1])
        end
    end
end

function stripe2cross!(elms::Vector{ApproxOperator.AbstractElement}, nds::Vector{𝑿ᵢ})
    nₑ = length(elms)
    nₚ = length(nds)
    nₛ = length(elms[1].𝓖)
    nᵢ = Int((nₑ/2)^0.5)
    data𝓒 = getfield(nds[1],:data)
    data𝓖 = getfield(elms[1].𝓖[1],:data)
    x = data𝓒[:x][2]
    y = data𝓒[:y][2]
    z = data𝓒[:z][2]
    push!(data𝓖[:x][2],zeros(nₑ*nₛ)...)
    push!(data𝓖[:y][2],zeros(nₑ*nₛ)...)
    push!(data𝓖[:z][2],zeros(nₑ*nₛ)...)
    push!(data𝓖[:𝑤][2],zeros(nₑ*nₛ)...)
    push!(data𝓖[:𝐴][2],zeros(nₑ)...)
    for i in 1:nᵢ
        for j in 1:nᵢ
            ind = 2*nᵢ*(j-1)+2*i-1
            x₁ = elms[ind].𝓒[1]
            x₂ = elms[ind].𝓒[2]
            x₃ = elms[ind+1].𝓒[3]
            x₄ = elms[ind+1].𝓒[1]
            nₚ += 1
            x₅ = 𝑿ᵢ((𝐼=nₚ,),data𝓒)
            push!(nds,x₅)
            push!(x,0.25*(x₁.x+x₂.x+x₃.x+x₄.x))
            push!(y,0.25*(x₁.y+x₂.y+x₃.y+x₄.y))
            push!(z,0.25*(x₁.z+x₂.z+x₃.z+x₄.z))
            elms[ind].𝓒[1:3] .= [x₁, x₂, x₅]
            resetx!(elms[ind])
            reset𝐴!(elms[ind])
            elms[ind+1].𝓒[1:3] .= [x₂, x₃, x₅]
            resetx!(elms[ind+1])
            reset𝐴!(elms[ind+1])

            𝓒 = [x₃, x₄, x₅]
            𝓖 = [𝑿ₛ((𝑔=k,𝐺=nₑ*nₛ+k,𝐶=nₑ+1,𝑠=3*(nₑ*nₛ+k)),data𝓖) for k in 1:nₛ]
            elm = Element{:Tri3}(𝓒,𝓖)
            resetx!(elm)
            reset𝐴!(elm)
            push!(elms,elm)
            nₑ+=1

            𝓒 = [x₄, x₁, x₅]
            𝓖 = [𝑿ₛ((𝑔=k,𝐺=nₑ*nₛ+k,𝐶=nₑ+1,𝑠=3*(nₑ*nₛ+k)),data𝓖) for k in 1:nₛ]
            elm = Element{:Tri3}(𝓒,𝓖)
            resetx!(elm)
            reset𝐴!(elm)
            push!(elms,elm)
            nₑ+=1
        end
    end
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

prescribe_Ωᵍᵖ = quote
    prescribe!(elements["Ωᵍᵖ"],:u=>(x,y,z)->u(x,y))
    prescribe!(elements["Ωᵍᵖ"],:v=>(x,y,z)->v(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
    prescribe!(elements["Ωᵍᵖ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
end