
using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf
include("input.jl")
# for i in 2:10
ndiv= 30
# ndiv_p=9
i=4


include("import_prescrible_ops.jl")
include("import_cook_membrane.jl")
include("vtk.jl")
elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes ,nodes_p = import_cook_membrane_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p = import_cantilever_T6P3("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_Q4P1("./msh/cantilever_quad_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑₚ = length(Ω)


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
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])




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
   
k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;fp]
d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
q  = d[2*nᵤ+1:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>q)
eval(VTK_mix_pressure)



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