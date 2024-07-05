
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv =3
nₚ = 5
# println(nₚ)
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")

nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωᵖ"])
set∇𝝭!(elements["Γ"])
D=1   #thermal conductivity coefficient
t=1 #thickness

n = 3
T(x,y) = (x+y)^n
∂T∂x(x,y) = n*(x+y)^abs(n-1)
∂T∂y(x,y) = n*(x+y)^abs(n-1)
∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
s(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))

eval(prescribe)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
       Operator{:∫∫Tᵢsᵢdxdy}(:t=>t),
       Operator{:T_error}(:D=>D,:t=>t),
]
opsᵛ = [
    Operator{:∫∫qᵢ∇Tⱼdxdy}(:t=>t),
]
opsᵈ = [
    Operator{:∫∫qᵢD⁻¹qⱼdxdy}(:D=>D,:t=>t),
]

kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(nₚ)


opsᵈ[1](elements["Ω"],kᵤᵤ)
opsᵛ[1](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
ops[3](elements["Ωᵖ"],f)
kₚᵤ⁻=kₚᵤ'*inv(kₚᵤ*kₚᵤ')
kₐ=-kᵤᵤ*kₚᵤ⁻
k=(kₐ'*inv(kₐ*kₐ'))*(kₚᵤ')
ops[2](elements["Γ"],kᵅ,f)
 q = (k+kᵅ)\f #temperatures

# k = [kᵤᵤ (kₚᵤ+kᵅ)';kₚᵤ+kᵅ kₚₚ]
# k = [kᵤᵤ kₚᵤ';kₚᵤ kₚₚ+kᵅ]
# f = [zeros(2*nᵤ);f+fᵅ]
# d = k\f
# d₁ = d[1:2:2*nᵤ] ##heat flux
# d₂ = d[2:2:2*nᵤ]
# q  = d[2*nᵤ+1:end]

# push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:T=>q)


set𝝭!(elements["Ωᵍᵖ"])
l2= ops[4](elements["Ωᵍᵖ"])
# h1,l2,h1_dil,h1_dev= ops[5](elements["Ωᵍ"],elements["Ωᵍᵖ"])
L2 = log10(l2)
# H1 = log10(h1)
# H1_dil = log10(h1_dil)
# H1_dev = log10(h1_dev)
           
println(L2)

# eval(VTK_mix_pressure)


