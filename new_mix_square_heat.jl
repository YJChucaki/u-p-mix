using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie
ndiv = 5
# ndiv_u = 11
nᵤ = 20

include("import_heat_conduction_infsup.jl")
include("wirteVTK.jl")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_quad_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_quad8_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
elements, nodes, nodes_u= import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nᵤ)*".msh")
# elements, nodes = import_cantilever_Q4P1("./msh/square_quad_"*string(ndiv)*".msh")
# elements, nodes = import_cantilever_Q8P3("./msh/square_quad8_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p = import_cantilever_mix("./msh/square_tri6_"*string(ndiv)*".msh","./msh/square_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p = import_cantilever_T6P3("./msh/square_tri6_"*string(ndiv)*".msh","./msh/square_"*string(ndiv)*".msh")   

nₚ = length(nodes)
nᵤ = length(nodes_u)
t = 1
D = 1
set𝝭!(elements["Ωᵘ"])
set𝝭!(elements["∂Ωᵘ"])
# set𝝭!(elements["Γ²ᵘ"])
set𝝭!(elements["Γ⁴ᵘ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["∂Ωᵖ"])
# set𝝭!(elements["Γ²ᵖ"])
set𝝭!(elements["Γ⁴ᵖ"])
##for Q4P1 
# nₚ = length(elements["Ωᵖ"])
##for Q8P3 
# nₚ = 3*length(elements["Ωᵖ"])
   


eval(prescribe)
prescribe!(elements["∂Ωᵘ"],:g=>(x,y,z)->0.0)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12,:t=>t),
       Operator{:∫vbdΩ}(),
       Operator{:L₂}(),
       Operator{:L₂_heat_flux}(),
]
opsᵛ = [
    Operator{:∫∫∇𝒑udxdy}(),
    Operator{:∫pᵢnᵢuds}(),
    Operator{:∫pᵢnᵢgⱼds}(),
]
opsᵈ = [
    Operator{:∫∫pᵢD⁻¹pⱼdxdy}(),
]

kₚₚ = zeros(2*nₚ,2*nₚ)
kₚₙ = zeros(2*nₚ,nᵤ)
kₚᵤ = zeros(2*nₚ,nᵤ)
kᵤᵤ = zeros(nᵤ,nᵤ)
fᵤ = zeros(nᵤ)
fₚ = zeros(2*nₚ)
 
opsᵈ[1](elements["Ωᵖ"],kₚₚ)
opsᵛ[1](elements["Ωᵖ"],elements["Ωᵘ"],kₚᵤ)
opsᵛ[2](elements["∂Ωᵖ"],elements["∂Ωᵘ"],kₚᵤ)
opsᵛ[3](elements["Γ⁴ᵖ"],elements["Γ⁴ᵘ"],kₚₙ,fₚ)
# ops[3](elements["Ωᵘ"],fᵤ)
kᵈ = -kₚᵤ'*(kₚₚ\kₚᵤ)
kᵛ = zeros(nᵤ,nᵤ)
vᵈ = eigvals(kᵈ)
vᵛ = eigvals(kᵛ)
v = eigvals(kᵛ,kᵈ)
γ = eigvals(kᵛ,kᵈ)
γ = eigvals(kᵈ,kᵛ)
# println(γ[2*nᵤ-nₚ+1])