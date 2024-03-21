using ApproxOperator, GLMakie, CairoMakie, TimerOutputs, JLD
import Gmsh: gmsh
ndiv= 8
 ndiv_p= 8
i=200

P = 1000
Ē = 3e6
ν̄ = 0.4999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 48
D = 12
I = D^3/12
EI = E*I

gmsh.initialize()
gmsh.open("./msh/cantilever_bubble_"*string(i)*".msh")
nodes_p = get𝑿ᵢ()
xᵖ = nodes_p.x
yᵖ = nodes_p.y
zᵖ = nodes_p.z
sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 1,γ = 2)
nₚ = length(nodes_p)
s =1.5*12/ndiv_p*ones(nₚ)
# s =3.5*0.5*ones(nₚ)
push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)
type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
q = Dict(load("jld/cantilever_mix_tri3_bubble_G3_"*string(i)*".jld"))
# q = Dict(load("jld/cantilever_mix_quad4_bubble_G3_"*string(i)*".jld"))
push!(nodes_p,:q=>q["q"])

gmsh.open("./msh/cantilever_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
integrationOrder = 1
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"],   integrationOrder, normal = true)
elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],   integrationOrder, normal = true)
elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],   integrationOrder, normal = true)

# nᵤ = length(nodes)


𝗠 = zeros(21)
ind = 20
xs = zeros(ind)
ys = zeros(ind)
color = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0, L, ind))
    for (J,ξ²) in enumerate(LinRange(-6.0, D/2, ind))
        indices = sp(ξ¹,ξ²,0.0)
        Nᵖ = zeros(length(indices))
        data = Dict([:x=>(1,[ξ¹]),:y=>(1,[ξ²]),:z=>(1,[0.0]),:𝝭=>(4,Nᵖ),:𝗠=>(0,𝗠)])
        𝓒 = [nodes_p[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
         p= 0.0
         for (i,xᵢ) in enumerate(𝓒)
            p  += Nᵖ[i]*xᵢ.q
        end 
        xs[I] = ξ¹
        ys[J] = ξ² 
        color[I,J] = p
    end
end
fig = Figure()

ax = Axis(fig[1, 1], aspect = 4)

hidespines!(ax)
hidedecorations!(ax)

# s=surface!(xs,ys, color, colormap=:coolwarm)
s = contourf!(xs,ys, color, colormap=:coolwarm)
Colorbar(fig[1, 2],s)

# elements
lwb = 2.5;lwm =2.5;mso =5;msx =15;ppu = 2.5;α = 0.7;
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
    lines!(x,y, linewidth = 0.3, color = :black)
end
# scatter!(x,y,marker = :circle, markersize = mso, color = :black)
lines!([0.0,L,L,0.0,0.0],[-D/2,-D/2,D/2,D/2,-D/2], linewidth = lwb, color = :black)
save("./png/cantilever_tri3_G3_"*string(i)*".png",fig)
# save("./png/cantilever_quad4_G3_"*string(i)*".png",fig)
# save("./png/cantilever_nomesh_"*string(i)*".png",fig)
fig