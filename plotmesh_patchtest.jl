using ApproxOperator,CairoMakie

include("import_patchtest.jl")


lwb = 2.5;lwm =2.5;mso =15;msx =15;ppu = 2.5;α = 0.7;
filename1 = "./msh/patchtest.msh"
# filename2 = "./msh/cantilever_bubble_72.msh"
savename = "./png/patchtest.png"
elms, nodes = import_patchtest("./msh/patchtest.msh")
# elms_p,~ = ApproxOperator.importmsh(filename2)

x = getfield(nodes[1],:data)[:x][2]
y = getfield(nodes[1],:data)[:y][2]
z = getfield(nodes[1],:data)[:z][2]
# xᵖ = elms_p["Ω"][1].x
# yᵖ = elms_p["Ω"][1].y
# zᵖ = elms_p["Ω"][1].z

if occursin("quad",filename1)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end

f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
hidespines!(ax)
hidedecorations!(ax)
lines!([0.0, 1.0, 1.0, 0.0, 0.0],[0.0, 0.0, 1.0, 1.0, 0.0], linewidth = lwb, color = :black)

for elm in elms["Ω"]
    id = [node.𝐼 for node in elm.𝓒]
    lines!(x[id[index]],y[id[index]], linewidth = lwm, color = :black)
end
scatter!(x,y,marker = :circle, markersize = mso, color = :black)

# for elm in elms_p["Ω"]
#     id = [i for i in elm.i]
#     # lines!(xᵖ[id[[1,2,3,1]]],yᵖ[id[[1,2,3,1]]], linewidth = lwm, color = :blue)
# end
# scatter!(xᵖ,yᵖ,marker = :xcross, markersize = msx, color = (:blue, α))
save(savename,f,px_per_unit = ppu)
f