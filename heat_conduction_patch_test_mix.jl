
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 5
ndiv_u = 3
# nᵤ = 5

# elements, nodes, nodes_u = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nᵤ)*".msh")
elements, nodes, nodes_u = import_patchtest_mix("./msh/patchtest_quad_"*string(ndiv)*".msh","./msh/patchtest_quad_"*string(ndiv_u)*".msh")
nₚ = length(nodes)
nᵤ = length(nodes_u)


set𝝭!(elements["Ωᵘ"])
set𝝭!(elements["∂Ωᵘ"])
set𝝭!(elements["Γᵘ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["∂Ωᵖ"])
set𝝭!(elements["Γᵖ"])
D=1   #thermal conductivity coefficient
t=1 #thickness

n =1
T(x,y) = (x+y)^n
∂T∂x(x,y) = n*(x+y)^abs(n-1)
∂T∂y(x,y) = n*(x+y)^abs(n-1)
P₁(x,y) = -∂T∂x(x,y)
P₂(x,y) = -∂T∂y(x,y)
∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           *(∂²T∂x²(x,y)+∂²T∂y²(x,y))

eval(prescribe)
prescribe!(elements["∂Ωᵘ"],:g=>(x,y,z)->0.0)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
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
    Operator{:∫∫pᵢD⁻¹pⱼdxdy}(:D=>D,:t=>t),
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
# opsᵛ[2](elements["Γᵖ"],elements["Γᵘ"],kₚᵤ)
opsᵛ[3](elements["Γᵖ"],elements["Γᵘ"],kₚₙ,fₚ)
ops[3](elements["Ωᵘ"],fᵤ)

# kₚᵤ = kₚᵤ+kₚₙ
# kᵈ = kₚₚ
# kᵛ = kₚᵤ*(kᵤᵤ\(-kₚᵤ'))
# vᵈ = eigvals(kᵈ)
# vᵛ = eigvals(kᵛ)
# v = eigvals(kᵛ,kᵈ)
# γ = eigvals(kᵛ,kᵈ)
# println(γ[2*nᵤ-nₚ+1])


k = [kₚₚ -kₚᵤ-kₚₙ;-kₚᵤ'-kₚₙ' kᵤᵤ]

# k = [kₚₚ -kₚᵤ;-kₚᵤ' kᵤᵤ]
f = [-fₚ;-fᵤ]
d = k\f
p₁ = d[1:2:2*nₚ] 
p₂ = d[2:2:2*nₚ]
u  = d[2*nₚ+1:end]

push!(nodes,:d₁=>p₁,:d₂=>p₂)
push!(nodes_u,:d=>u)


set𝝭!(elements["Ωᵍᵘ"])
set𝝭!(elements["Ωᵍᵖ"])
l2_u= ops[4](elements["Ωᵍᵘ"])
l2_p= ops[5](elements["Ωᵍᵖ"])
L2_u = log10(l2_u)
L2_p = log10(l2_p)

           
println(L2_u)
println(L2_p)

# eval(VTK_mix_pressure)

dₚᵤ = zeros(2*nₚ + nᵤ)
dₚ = zeros(2*nₚ)
dᵤ = zeros(nᵤ)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    dₚᵤ[2*i-1] = -∂T∂x(x,y)
    dₚᵤ[2*i]   = -∂T∂y(x,y)
    dₚ[2*i-1] = -∂T∂x(x,y)
    dₚ[2*i]   = -∂T∂y(x,y)
end
for (i,node) in enumerate(nodes_u)
    x = node.x
    y = node.y
    dₚᵤ[2*nₚ+i] = T(x,y)
    dᵤ[i] = T(x,y)
end

err1 = kₚₚ*dₚ - kₚᵤ*dᵤ
err2 = kₚₙ*dᵤ - fₚ
err3 = (kₚᵤ+kₚₙ)'*dₚ-fᵤ
err4 = kₚₚ*dₚ - (kₚᵤ+kₚₙ)*dᵤ + fₚ
err5 = k*dₚᵤ-f
err6 = k*d-f
