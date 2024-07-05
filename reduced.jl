using  ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf,Pardiso
include("import_prescrible_ops.jl")                       
include("import_cantilever.jl")
include("wirteVTK.jl")

ndiv= 8
ndiv_p=8
# elements,nodes,nodes_p= import_quad_GI1("./msh/square_quad_"*string(ndiv)*".msh","./msh/square_quad_"*string(ndiv_p)*".msh")
# elements,nodes,nodes_p= import_quad_GI1("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
# elements,nodes,nodes_p= import_fem_tri3_GI1("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
elements,nodes,nodes_p= import_cantilever_reduce("./msh/cantilever_quad_"*string(ndiv)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)
P = 1000
Ē = 3e6
ν̄ =0.49999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
L = 48
D = 12
I = D^3/12
EI = E*I

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵛ"])
set∇𝝭!(elements["Ωᵛ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
eval(prescribeForGauss)
eval(prescribeForPenalty)


ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e12*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:Hₑ_Incompressible}(:E=>Ē,:ν=>ν̄ ),

]
opsᵛ = [
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]

kᵛ = zeros(2*nᵤ,2*nᵤ)
kᵈ = zeros(2*nᵤ,2*nᵤ)
kᵍ = zeros(2*nₚ,2*nₚ)
f = zeros(2*nᵤ)

opsᵛ[1](elements["Ωᵛ"],kᵛ)
opsᵈ[1](elements["Ω"],kᵈ)  
# ops[2](elements["Ω"],f)
ops[3](elements["Γᵗ"],f)
ops[4](elements["Γᵍ"],kᵍ,f)

d = (kᵛ+kᵈ+kᵍ)\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ωᵍ"])
h1,l2 = ops[6](elements["Ωᵍ"])
L2 = log10(l2)
H1 = log10(h1)
println(L2,H1)
