using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie
# NP=[40,80,120,140]
# for n=1:4
    # i=NP[n]
ndiv= 8
#  ndiv_p= 8
i=260
# 40,60-3
# 80-4
# 100,120-5
# 160,200-7

include("import_prescrible_ops.jl")
include("import_cantilever.jl")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
# elements, nodes ,nodes_p = import_cantilever_mix_quad4("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix_tri3("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix_quad4("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
    nᵤ = length(nodes)
    nₚ = length(nodes_p)
    nₘ=21
    P = 1000
    Ē = 3e6
    # Ē = 1.0
    ν̄ = 0.4999999
    # ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 48
    D = 12
    I = D^3/12
    EI = E*I
    K=Ē/3/(1-2ν̄ )
    eval(prescribeForGauss)
    eval(prescribeForPenalty)

    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍᵖ"])
   

    

    eval(opsupmix)
    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    kᵤₚ = zeros(2*nᵤ,nₚ)
    kₚₚ = zeros(nₚ,nₚ)
    f = zeros(2*nᵤ)
    fp= zeros(nₚ)
    opsup[3](elements["Ω"],kᵤᵤ)
    opsup[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
    opsup[5](elements["Ωᵖ"],kₚₚ)
    opsup[6](elements["Γᵗ"],f)
    αᵥ = 1e9

    eval(opsPenalty)
    opsα[1](elements["Γᵍ"],kᵤᵤ,f)
    opsα[2](elements["Γᵍ"],elements["Γᵍᵖ"],kᵤₚ,fp)

    #  kₚₚ⁻¹=inv(kₚₚ)
    # d = (kᵤᵤ-kᵤₚ*kₚₚ⁻¹*kᵤₚ')\f
    # q=-kₚₚ⁻¹*kᵤₚ'*d
    # d₃ = d[1:2*nᵤ]
    # d = (kᵤᵤ-kᵤₚ*kₚₚ⁻¹*kᵤₚ')\f
    # q=-kₚₚ⁻¹*kᵤₚ'*d

    k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
    f = [f;fp]
    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    q  = d[2*nᵤ+1:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)

    h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = opsup[8](elements["Ω"],elements["Ωᵖ"])
    L2 = log10(l2)
    H1 = log10(h1)
    H1_dil = log10(h1_dil)
    H1_dev = log10(h1_dev)
   
    println(L2,H1)
    println(H1_dil,H1_dev)
    println(l2,h1)
    println(h1_dil,h1_dev)
    # h = log10(10.0/ndiv)

#     index = 40:50
#     XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
#         Sheet = xf[2]
#         ind = findfirst(n->n==ndiv,index)+1
#         Sheet["F"*string(ind)] = h
#         Sheet["G"*string(ind)] = L2
#         Sheet["H"*string(ind)] = H1

# @save compress=true "jld/cantilever_mix_tri3_"*string(ndiv)*".jld" q
# @save compress=true "jld/cantilever_mix_tri3_bubble_G30_"*string(i)*".jld" q
# @save compress=true "jld/cantilever_mix_quad4_bubble_G3_"*string(i)*".jld" q
# @save compress=true "jld/cantilever_mix_quad4_"*string(ndiv)*".jld" q
#     end
# end

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
# s = contourf!(xs,ys, color, colormap=:coolwarm,levels=-1000:200:1000)
s = contourf!(xs,ys, color, colormap=:coolwarm)
Colorbar(fig[1, 2], s)

# # elements
lwb = 2.5;lwm =2.5;mso =5;msx =15;ppu = 2.5;α = 0.7;
for elm in elements["Ω"]
   
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
   
    lines!(x,y, linewidth = 0.3, color = :black)

end
# scatter!(x,y,marker = :circle, markersize = mso, color = :black)
lines!([0.0,L,L,0.0,0.0],[-D/2,-D/2,D/2,D/2,-D/2], linewidth = lwb, color = :black)
# save("./png/cantilever_"*string(i)*".png",fig)
# save("./png/cantilever_tri3_G3_level_"*string(i)*".png",fig)
# save("./png/cantilever_tri3_G3_nonunoform_level_"*string(i)*".png",fig)
fig
# end