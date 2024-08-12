using ApproxOperator,CairoMakie,Tensors, BenchmarkExample, Statistics
import Gmsh: gmsh
lwb = 1.5;lwm =0.1;mso =3;msx =2;ppu = 2.5;α = 0.7;
# filename1 = "./msh/cantilever_32.msh"
# filename2 = "./msh/cantilever_bubble_4165.msh"
# savename = "./png/cantilever_32_4165.png"
filename1 = "./msh/plate_with_hole_45.msh"
filename2 = "./msh/plate_with_hole_2.msh"
savename = "./png/plate_with_hole_2_45.png"

gmsh.initialize()
gmsh.open(filename1)

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
integrationOrder = 2
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)


gmsh.open(filename2)
nodes_p = get𝑿ᵢ()
xᵖ = nodes_p.x
yᵖ = nodes_p.y
zᵖ = nodes_p.z
elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"])


if occursin("quad",filename1)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end

f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
hidespines!(ax)
hidedecorations!(ax)
L = 48.
b = 12.

for elm in elements["Ωᵖ"]
    id = [node.𝐼 for node in elm.𝓒]
    lines!(xᵖ[id[index]],yᵖ[id[index]], linewidth = lwm, color = :grey)
end
# lines!([0.0,L,L,0.0,0.0],[-b/2,-b/2,b/2,b/2,-b/2], linewidth = lwb, color = :black)
scatter!(xᵖ,yᵖ,marker = :circle, markersize = mso, color = :black)

for elm in elements["Ωᵖ"]
    id = [node.𝐼 for node in elm.𝓒]
    # lines!(xᵖ[id[[1,2,3,1]]],yᵖ[id[[1,2,3,1]]], linewidth = lwm, color = :blue)
end
scatter!(x,y,marker = :xcross, markersize = msx, color = (:blue, α))

save(savename,f,px_per_unit = ppu)
f