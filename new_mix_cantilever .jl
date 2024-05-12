using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie

ndiv= 33
i=200

include("import_prescrible_ops.jl")                       
include("import_cantilever.jl")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p = import_cantilever_T6P3("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_Q4P1("./msh/cantilever_quad_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_Q4R1("./msh/cantilever_quad_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
elements, nodes  = import_cantilever_Q8P3("./msh/cantilever_quad8_"*string(ndiv)*".msh")
    nᵤ = length(nodes)
    nₚ = length(nodes_p)
    nₑᵤ = length(elements["Ω"])
    ##for Q4P1 
    # nₚ = length(elements["Ωᵖ"])
    ##for Q8P3
    # nₚ = 3*length(elements["Ωᵖ"])
    P = 1000
    Ē = 3e6
    # Ē = 1.0
    # ν̄ = 0.4999999
    ν̄ = 0.3
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
    # set𝝭!(elements["Γᵍᵖ"])
   

    

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
    # opsα[2](elements["Γᵍ"],elements["Γᵍᵖ"],kᵤₚ,fp)

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
    # push!(nodes_p,:q=>q)

    # h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = opsup[8](elements["Ωᵍ"],elements["Ωᵖ"])
    h1,l2 = opsup[9](elements["Ωᵍ"])
    L2 = log10(l2)
    H1 = log10(h1)
    # H1_dil = log10(h1_dil)
    # H1_dev = log10(h1_dev)
   
    println(L2,H1)
    # println(H1_dil,H1_dev)
    # println(l2,h1)
    # println(h1_dil,h1_dev)
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

##Saved as VTK


# fo = open("./vtk/cook_membrance_rkgsi_mix_"*string(ndiv_𝑢)*".vtk","w")
fo = open("./vtk/cantilever_tri3_mix_"*string(ndiv_𝑢)*".vtk","w")
@printf fo "# vtk DataFile Version 2.0\n"
@printf fo "cantilever_tri3_mix\n"
@printf fo "ASCII\n"
@printf fo "DATASET POLYDATA\n"
@printf fo "POINTS %i float\n" nᵤ
for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
end
@printf fo "POLYGONS %i %i\n" nₑᵤ 4*nₑᵤ
for ap in elements["Ω"]
    𝓒 = ap.𝓒s
    @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
end
@printf fo "POINT_DATA %i\n" nᵤ
@printf fo "VECTORS U float\n"
for p in nodes
   
    @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
end

@printf fo "TENSORS STRESS float\n"
for p in elements["Ω"]
    𝓒 = p.𝓒
    𝓖 = p.𝓖
    ε₁₁ = 0.0
    ε₂₂ = 0.0
    ε₁₂ = 0.0

    for (i,ξ) in enumerate(𝓖)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (j,xⱼ) in enumerate(𝓒)
            ε₁₁ += B₁[j]*xⱼ.d₁
            ε₂₂ += B₂[j]*xⱼ.d₂
            ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
        end
    end
    σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
    σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
    σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
    @printf fo "%f %f %f\n" σ₁₁ σ₁₂ 0.0
    @printf fo "%f %f %f\n" σ₁₂ σ₂₂ 0.0
    @printf fo "%f %f %f\n" 0.0 0.0 0.0
end
close(fo)

 ##contour!
# 𝗠 = zeros(21)
# ind = 20
# xs = zeros(ind)
# ys = zeros(ind)
# color = zeros(ind,ind)

# for (I,ξ¹) in enumerate(LinRange(0.0, L, ind))
#     for (J,ξ²) in enumerate(LinRange(-6.0, D/2, ind))
#         indices = sp(ξ¹,ξ²,0.0)
#         Nᵖ = zeros(length(indices))
#         data = Dict([:x=>(1,[ξ¹]),:y=>(1,[ξ²]),:z=>(1,[0.0]),:𝝭=>(4,Nᵖ),:𝗠=>(0,𝗠)])
#         𝓒 = [nodes_p[k] for k in indices]
#         𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
#         ap = type(𝓒,𝓖)
#         set𝝭!(ap)
#          p= 0.0       
#         for (i,xᵢ) in enumerate(𝓒)
#             p  += Nᵖ[i]*xᵢ.q
           
#         end 
#         xs[I] = ξ¹
#         ys[J] = ξ² 
#         color[I,J] = p
        
#     end
# end

# fig = Figure()
# ax = Axis(fig[1, 1], aspect = 4)
# hidespines!(ax)
# hidedecorations!(ax)

# # s=surface!(xs,ys, color, colormap=:coolwarm)
# # s = contourf!(xs,ys, color, colormap=:coolwarm,levels=-1000:200:1000)
# s = contourf!(xs,ys, color, colormap=:coolwarm)
# Colorbar(fig[1, 2], s)

# # # elements
# lwb = 2.5;lwm =2.5;mso =5;msx =15;ppu = 2.5;α = 0.7;
# for elm in elements["Ω"]
   
#     x = [x.x for x in elm.𝓒[[1,2,3,1]]]
#     y = [x.y for x in elm.𝓒[[1,2,3,1]]]
   
#     lines!(x,y, linewidth = 0.3, color = :black)

# end
# # scatter!(x,y,marker = :circle, markersize = mso, color = :black)
# lines!([0.0,L,L,0.0,0.0],[-D/2,-D/2,D/2,D/2,-D/2], linewidth = lwb, color = :black)
# # save("./png/cantilever_"*string(i)*".png",fig)
# # save("./png/cantilever_tri3_G3_level_"*string(i)*".png",fig)
# # save("./png/cantilever_tri3_G3_nonunoform_level_"*string(i)*".png",fig)
# # save("./png/cantilever_tri6_G3_level_"*string(i)*".png",fig)
# fig
