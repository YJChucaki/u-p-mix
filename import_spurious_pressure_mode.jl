
import Gmsh: gmsh
using Statistics

function import_test(filename)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    integrationOrder = 8
    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    return elements,nodes
end

function import_test_2(filename1::String, filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    x = nodes_p.x
    y = nodes_p.y
    z = nodes_p.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)

    Ω = getElements(nodes_p, entities["Ω"])
    s, ~ = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    integrationOrder = 8
    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    elements["Ωₚ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder, sp)
    𝗠 = (0,zeros(21))
    push!(elements["Ωₚ"], :𝝭=>:𝑠)
    push!(elements["Ωₚ"], :𝗠=>𝗠)

    return elements,nodes,nodes_p
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

