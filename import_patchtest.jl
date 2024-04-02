using Gmsh
using CairoMakie
# using GLMakie

function import_patchtest_stripe(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    integrationOrder = 2
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)

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

    integrationOrder = 2
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    stripe2cross!(elements["Ω"], nodes)

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