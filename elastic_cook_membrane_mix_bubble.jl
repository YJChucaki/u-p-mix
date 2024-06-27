
using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf
include("input.jl")
# for i in 2:10
ndiv= 15
# ndiv_p=9
i=10


include("import_prescrible_ops.jl")
include("import_cook_membrane.jl")
include("wirteVTK.jl")
elements, nodes ,nodes_p,Ω = import_cook_membrane_mix_bubble("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_quad_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes = import_cook_membrane_Q4P1("./msh/cook_membrane_quad_"*string(ndiv)*".msh")
nᵤ = length(nodes)
nₒ = length(elements["Ω"])
# nₚ = length(nodes_p)
# nₑₚ = length(Ω)
nₚ = length(elements["Ωᵖ"])
nₑ = length(elements["Ω"])
# κ = 400942
# μ = 80.1938
# E = 9*κ*μ/(3*κ+μ)
# ν = (3*κ-2*μ)/2/(3*κ+μ)
Ē = 70.0
# ν = 0.3333
ν̄  =0.4999999
E = Ē/(1.0-ν̄^2)
 ν = ν̄/(1.0-ν̄)
eval(prescribeForPenalty)
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
    # fₛ = zeros(4*nₛ)
    fₒ = zeros(2*nₒ)
    dᵤ = zeros(2*nᵤ)
    dₚ = zeros(nₚ)
    # dₛ = zeros(4*nₛ)
        
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
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)
    
# push!(nodes_p,:q=>q)
# eval(VTK_mix_pressure)
# eval(VTK_Q4P1_displacement_pressure)


# # fo = open("./vtk/cook_membrance_rkgsi_mix_"*string(ndiv_𝑢)*".vtk","w")
# fo = open("./vtk/cook_membrance_rkgsi_"*string(ndiv_𝑢)*".vtk","w")
# @printf fo "# vtk DataFile Version 2.0\n"
# @printf fo "cook_membrance_rkgsi_mix\n"
# @printf fo "ASCII\n"
# @printf fo "DATASET POLYDATA\n"
# @printf fo "POINTS %i float\n" nₚ
# for p in nodes
#     @printf fo "%f %f %f\n" p.x p.y p.z
# end
# @printf fo "POLYGONS %i %i\n" nₑ 4*nₑ
# for ap in elms["Ω"]
#     𝓒 = ap.vertices
#     @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
# end
# @printf fo "POINT_DATA %i\n" nₚ
# @printf fo "VECTORS U float\n"
# for p in elements["Ωᶜ"]
#     ξ = collect(p.𝓖)[1]
#     N = ξ[:𝝭]
#     u₁ = 0.0
#     u₂ = 0.0
#     for (i,x) in enumerate(p.𝓒)
#         u₁ += N[i]*x.d₁
#         u₂ += N[i]*x.d₂
#     end
#     @printf fo "%f %f %f\n" u₁ u₂ 0.0
# end

# @printf fo "TENSORS STRESS float\n"
# for p in elements["Ωᶜ"]
#     𝓒 = p.𝓒
#     𝓖 = p.𝓖
#     ε₁₁ = 0.0
#     ε₂₂ = 0.0
#     ε₁₂ = 0.0

#     for (i,ξ) in enumerate(𝓖)
#         B₁ = ξ[:∂𝝭∂x]
#         B₂ = ξ[:∂𝝭∂y]
#         for (j,xⱼ) in enumerate(𝓒)
#             ε₁₁ += B₁[j]*xⱼ.d₁
#             ε₂₂ += B₂[j]*xⱼ.d₂
#             ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
#         end
#     end
#     σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
#     σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
#     σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
#     @printf fo "%f %f %f\n" σ₁₁ σ₁₂ 0.0
#     @printf fo "%f %f %f\n" σ₁₂ σ₂₂ 0.0
#     @printf fo "%f %f %f\n" 0.0 0.0 0.0
# end
# close(fo)

a = elements["Ω"][end]
ξs = collect(a.𝓖)
𝝭 = ξs[3][:𝝭]
u₂ = 0.0
for (i,x) in enumerate(a.𝓒)
    global u₂ += 𝝭[i]*x.d₂
end
h = nᵤ/nₚ

println(u₂)
# println(nₚ)
# index = 10:30
#     XLSX.openxlsx("./xlsx/cook.xlsx", mode="rw") do xf
#         Sheet = xf[5]
#         ind = findfirst(n->n==ndiv_p,index)+1
#         Sheet["B"*string(ind)] = h
#         Sheet["C"*string(ind)] = u₂
       
#     end
# end