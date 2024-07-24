using ApproxOperator, JLD,LinearAlgebra, Printf ,Pardiso, Tensors

ndiv= 33
i= 1042
# ndiv_p=4
include("import_prescrible_ops.jl")                       
include("import_cantilever.jl")
include("wirteVTK.jl")

# elements, nodes, Ω  = import_cantilever_mix_HR("./msh/cantilever.msh","./msh/cantilever_bubble_"*string(i)*".msh")
elements, nodes, nodes_p, Ω  = import_cantilever_mix_bubble("./msh/square_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes, nodes_p, Ω  = import_cantilever_mix_bubble("./msh/cantilever_HR_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes, nodes_p, Ω  = import_cantilever_mix_bubble("./msh/cantilever_HR_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")   
    nₒ = length(elements["Ω"])
    nₑ = length(elements["Ω"])
    nᵤ = length(nodes)
    nₚ = length(nodes_p)

    nₛ = 3*nₑ
    nₑₚ = length(Ω)
    ##for Q4P1 
    # nₚ = length(elements["Ωᵖ"])
    ##for Q8P3
    # nₚ = 3*length(elements["Ωᵖ"])
    P = 1000
    Ē = 3e6
    # Ē = 1.0
    # ν̄ = 0.499999999
    ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 1
    D = 1
    I = D^3/12
    EI = E*I
    K=Ē/3/(1-2ν̄ )
    # eval(prescribeForGauss)
    # eval(prescribeForPenalty)

    eval(prescribeForSquare)
    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    # set𝝭!(elements["Ωˢ"])
    set𝝭!(elements["Ωᵖ"])
    set∇𝝭!(elements["Ωᵇ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵖ"])
    # set𝝭!(elements["Γˢ"])
    # set𝝭!(elements["Γᵍᵖ"])
   
opsᵖ = [
    Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫pnᵢgᵢds}(),
]

opsˢ = [
    Operator{:∫∫δsᵢⱼsᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫∫sᵢⱼεᵢⱼdxdy}(),
    Operator{:∫sᵢⱼnⱼgᵢds}(),
]

ops = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢgᵢds}(:α=>1e10*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄),
]
  
    
kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
# kₛᵤ = zeros(4*nₛ,2*nᵤ)
kₒᵤ = zeros(2*nₒ,2*nᵤ)
kₚₒ = zeros(nₚ,2*nₒ)
kₚₚ = zeros(nₚ,nₚ)
# kₛₛ = zeros(4*nₛ,4*nₛ)
kₒₒ = zeros(2*nₒ,2*nₒ)
fᵤ = zeros(2*nᵤ)
fₚ = zeros(nₚ)
fₛ = zeros(4*nₛ)
fₒ = zeros(2*nₒ)
dᵤ = zeros(2*nᵤ)
dₚ = zeros(nₚ)
dₛ = zeros(4*nₛ)
    
opsᵖ[1](elements["Ωᵖ"],kₚₚ)
opsᵖ[2](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
opsᵖ[2](elements["Ωᵇ"],elements["Ωᵖ"],kₚₒ)
# opsᵖ[3](elements["Γᵍ"],elements["Γᵖ"],kₚᵤ,fₚ)

# opsˢ[1](elements["Ωˢ"],kₛₛ)
# opsˢ[2](elements["Ω"],elements["Ωˢ"],kₛᵤ)
# opsˢ[3](elements["Γᵍ"],elements["Γˢ"],kₛᵤ,fₛ)



ops[1](elements["Ω"],kᵤᵤ)
ops[1](elements["Ωᵇ"],kₒₒ)
ops[1](elements["Ω"],elements["Ωᵇ"],kₒᵤ)
ops[2](elements["Γᵗ"],fᵤ)
ops[4](elements["Γᵍ"],kᵤᵤ,fᵤ)



# k = [kᵤᵤ kₚᵤ' kₛᵤ' kₒᵤ';
#      kₚᵤ kₚₚ zeros(nₚ,4*nₛ) zeros(nₚ,2*nₒ);
#      kₛᵤ zeros(4*nₛ,nₚ) kₛₛ zeros(4*nₛ,2*nₒ);
#      kₒᵤ zeros(2*nₒ,nₚ) zeros(2*nₒ,4*nₛ) kₒₒ]
k = [kᵤᵤ kₚᵤ' kₒᵤ';
     kₚᵤ kₚₚ kₚₒ;
     kₒᵤ kₚₒ' kₒₒ]
# f = [fᵤ;fₚ;fₛ;fₒ]
f = [fᵤ;fₚ;fₒ]
    d = k\f
 

d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
q  = d[2*nᵤ+1:2*nᵤ+nₚ]
# push!(nodes,:d₁=>d₁,:d₂=>d₂)
# push!(nodes_p,:q=>q)
# u = d[1:2*nᵤ]
# d = [u;q]
# f = [fᵤ;fₚ]
# kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
# kᵤₚ = zeros(2*nᵤ,nₚ)
# kₚₚ = zeros(nₚ,nₚ)

# d⁻=d'*inv(d*d')
# kₑ = f*d⁻
# kₑ = f/d
# for I in 1:2*nᵤ
#     for J in 1:2*nᵤ
#         kᵤᵤ[I,J] =kₑ[I,J]
#     end
# end 

# for I in 1:2*nᵤ
#     for J in  2*nᵤ+1:2*nᵤ+nₚ
#         kᵤₚ[I,J-2*nᵤ] = kₑ[I,J]
#     end
# end 

# for I in 2*nᵤ+1:2*nᵤ+nₚ
#     for J in 2*nᵤ+1:2*nᵤ+nₚ
#         kₚₚ[I-2*nᵤ,J-2*nᵤ] = kₑ[I,J]
#     end
# end 

kᵈ = kᵤᵤ- kₒᵤ'*(kₒₒ\kₒᵤ)
kᵛ =-(kₚᵤ'+kₒᵤ'*(kₒₒ\kₚₒ'))*((kₚₚ+kₚₒ*(kₒₒ\kₚₒ'))\(kₚᵤ-kₚₒ*(kₒₒ\kₒᵤ)))
vᵈ = eigvals(kᵈ)
vᵛ = eigvals(kᵛ)
γ = eigvals(kᵛ,kᵈ)
println(γ[2*nᵤ-nₚ+1])   



    # h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = ops[6](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = ops[9](elements["Ωᵍ"])
    # L2 = log10(l2)
    # H1 = log10(h1)
    # H1_dil = log10(h1_dil)
    # H1_dev = log10(h1_dev)
   
    # println(L2,H1)
    # println(H1_dil,H1_dev)
    # println(l2,h1)
    # println(h1_dil,h1_dev)
    # h = log10(10.0/ndiv)

    
    # eval(VTK_mix_pressure)
    # eval(VTK_mix_pressure_u)
    # eval(VTK_mix_displacement)
    # eval(VTK_Q4P1_displacement_pressure)
    # eval(VTK_T6P3_pressure)



#  #contour!
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

# # # # elements
# lwb = 2.5;lwm =2.5;mso =5;msx =15;ppu = 2.5;α = 0.7;
# # for elm in elements["Ω"]
   
# #     x = [x.x for x in elm.𝓒[[1,2,3,1]]]
# #     y = [x.y for x in elm.𝓒[[1,2,3,1]]]
   
# #     lines!(x,y, linewidth = 0.3, color = :black)

# # end
# # scatter!(x,y,marker = :circle, markersize = mso, color = :black)
# lines!([0.0,L,L,0.0,0.0],[-D/2,-D/2,D/2,D/2,-D/2], linewidth = lwb, color = :black)
# # save("./png/cantilever_"*string(i)*".png",fig)
# # save("./png/cantilever_tri3_G3_level_"*string(i)*".png",fig)
# # save("./png/cantilever_tri3_G3_nonunoform_level_"*string(i)*".png",fig)
# # save("./png/cantilever_tri6_G3_level_"*string(i)*".png",fig)
# fig
