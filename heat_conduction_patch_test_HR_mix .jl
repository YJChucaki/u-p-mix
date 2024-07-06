
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv =5
nₚ = 15
# println(nₚ)
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")

nₚ = length(nodes)
nᵤ = length(nodes_p)
# nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set∇𝝭!(elements["Ωᵘ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵖ"])
set∇𝝭!(elements["Γᵘ"])
D=1   #thermal conductivity coefficient
t=1 #thickness

n =1
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
    Operator{:∫qᵢnᵢgⱼds}(:t=>t),
]
opsᵈ = [
    Operator{:∫∫qᵢD⁻¹qⱼdxdy}(:D=>D,:t=>t),
]

kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kᵤᵤ = zeros(2*nₚ,2*nₚ)
kₚₙ = zeros(2*nₚ,nᵤ)
kₚᵤ = zeros(nᵤ,2*nₚ)
kₚₚ = zeros(nᵤ,nᵤ)
f = zeros(nᵤ)
fₚ = zeros(2*nₚ)


opsᵈ[1](elements["Ωᵖ"],kᵤᵤ)
opsᵛ[1](elements["Ωᵖ"],elements["Ωᵘ"],kₚᵤ)
opsᵛ[2](elements["Γᵖ"],elements["Γᵘ"],kₚₙ,fₚ)

ops[3](elements["Ωᵘ"],f)
# kₚᵤ⁻=kₚᵤ'*inv(kₚᵤ*kₚᵤ')
# kₐ=-kᵤᵤ*kₚᵤ⁻
# k=(kₐ'*inv(kₐ*kₐ'))*(kₚᵤ')
# ops[2](elements["Γᵘ"],kᵅ,fᵅ)
#  q = (k+kᵅ)\(f+fᵅ) #temperatures


# k = [kᵤᵤ (kₚᵤ+kᵅ)';kₚᵤ+kᵅ kₚₚ]
k = [kᵤᵤ -kₚᵤ'+kₚₙ;-kₚᵤ-kₚₙ' kₚₚ]
f = [fₚ;f]
d = k\f
p₁ = d[1:2:2*nₚ] ##heat flux
p₂ = d[2:2:2*nₚ]
u  = d[2*nₚ+1:end]

push!(nodes,:d₁=>p₁,:d₂=>p₂)
push!(nodes_p,:T=>u)


set𝝭!(elements["Ωᵍᵘ"])
l2= ops[4](elements["Ωᵍᵘ"])
# h1,l2,h1_dil,h1_dev= ops[5](elements["Ωᵍ"],elements["Ωᵍᵖ"])
L2 = log10(l2)
# H1 = log10(h1)
# H1_dil = log10(h1_dil)
# H1_dev = log10(h1_dev)
           
println(L2)

# eval(VTK_mix_pressure)


