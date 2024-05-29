using ApproxOperator, JLD, CairoMakie 

import Gmsh: gmsh

include("import_patchtest.jl")

ndiv = 90
α = 0.0

## import nodes
savename = "./png/patchtest_bubble_"*string(ndiv)*"_c.png"
gmsh.initialize()
gmsh.open("msh/patchtest_bubble_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
Ω = getElements(nodes, entities["Ω"])
gmsh.finalize()
s, var𝐴 = cal_area_support(Ω)
s = 1.5*s*ones(length(nodes))
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
𝗠 = zeros(6)

# data = load("jld/stability_tri3_"*string(ndiv)*".jld")
# d₁ = value(data[:u])
d₁, d₂, p = load("jld/stability_tri3_"*string(ndiv)*".jld")
push!(nodes,:d₁=>d₁[2],:d₂=>d₂[2],:p=>p[2])

ind = 10
xs = zeros(ind)
ys = zeros(ind)
zs = zeros(ind,ind)
color = zeros(ind,ind)
for (I,xₛ) in enumerate(LinRange(0.0, 1.0, ind))
    for (J,yₛ) in enumerate(LinRange(0.0, 1.0, ind))
        indices = sp(xₛ,yₛ,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[xₛ]),:y=>(2,[yₛ]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        p = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            # u₁ += N[i]*xᵢ.d₁
            # u₂ += N[i]*xᵢ.d₂
            p += N[i]*xᵢ.p
        end
        xs[I] = xₛ + α*u₁
        ys[J] = yₛ + α*u₂
        zs[I,J] = 0.0
        color[I,J] = p
    end
end

fig = Figure()

ax = Axis(fig[1, 1],aspect = DataAspect())

hidespines!(ax)
hidedecorations!(ax)
# lines!([Point(0,0,0), Point(1,0,0), Point(1,1,0), Point(0,1,0), Point(0,0,0)],color=:black)
s = surface!(ax,xs,ys,zs, colorrange = (0, 1), color=color, colormap=:coolwarm)
# s = surface!(ax,xs,ys,zs, color=color, colormap=:coolwarm)
# s = surface!(ax,xs,ys,color)
# Colorbar(fig[2, 1], s, vertical = false)
Colorbar(fig[1, 2], s, vertical = true)
save(savename,fig,px_per_unit = 2.5)
fig
