
using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf
include("input.jl")
# for i in 2:10
ndiv= 25



include("import_prescrible_ops.jl")
include("import_cook_membrane.jl")
include("wirteVTK.jl")
elements, nodes  = import_cook_membrane_fem("./msh/cook_membrane_tri6_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_quad_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
nᵤ = length(nodes)
nₑ = length(elements["Ω"])


# κ = 400942
# μ = 80.1938
# E = 9*κ*μ/(3*κ+μ)
# ν = (3*κ-2*μ)/2/(3*κ+μ)
Ē = 70.0
# ν = 0.3333
ν̄  =0.4999999999
E = Ē/(1.0-ν̄^2)
 ν = ν̄/(1.0-ν̄)
eval(prescribeForPenalty)

set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])


ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e12*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
    
]
k = zeros(2*nᵤ,2*nᵤ)
f = zeros(2*nᵤ)

ops[1](elements["Ω"],k)
 ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],k,f)

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
# eval(VTK_displacement) 



a = elements["Ω"][end]
ξs = collect(a.𝓖)
𝝭 = ξs[3][:𝝭]
u₂ = 0.0
for (i,x) in enumerate(a.𝓒)
    global u₂ += 𝝭[i]*x.d₂
end
# h = nᵤ/nₚ

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