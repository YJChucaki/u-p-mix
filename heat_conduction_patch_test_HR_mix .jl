
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 3
ndiv2 = 11
# nₚ = 15
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv2)*".msh")

nₚ = length(nodes_p)
nᵤ = length(nodes)
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
b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))

eval(prescribe)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
       Operator{:∫vbdΩ}(),
       Operator{:L₂}(),
]
opsᵛ = [
    Operator{:∫∫pᵢ∇uⱼdxdy}(),
    Operator{:∫pᵢnᵢgⱼds}(),
]
opsᵈ = [
    Operator{:∫∫pᵢD⁻¹pⱼdxdy}(:D=>D,:t=>t),
]

kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kₚₚ = zeros(2*nₚ,2*nₚ)
kₚₙ = zeros(2*nₚ,nᵤ)
kₚᵤ = zeros(nᵤ,2*nₚ)
kᵤᵤ = zeros(nᵤ,nᵤ)
fᵤ = zeros(nᵤ)
fₚ = zeros(2*nₚ)


opsᵈ[1](elements["Ωᵖ"],kₚₚ)
opsᵛ[1](elements["Ωᵖ"],elements["Ωᵘ"],kₚᵤ)
opsᵛ[2](elements["Γᵖ"],elements["Γᵘ"],kₚₙ,fₚ)
ops[3](elements["Ωᵘ"],fᵤ)




k = [kₚₚ -kₚᵤ'-kₚₙ;-kₚᵤ-kₚₙ' kᵤᵤ]
f = [-fₚ;fᵤ]
d = k\f
p₁ = d[1:2:2*nₚ] 
p₂ = d[2:2:2*nₚ]
u  = d[2*nₚ+1:end]

push!(nodes_p,:d₁=>p₁,:d₂=>p₂)
push!(nodes,:d=>u)


set𝝭!(elements["Ωᵍᵘ"])
l2= ops[4](elements["Ωᵍᵘ"])
L2 = log10(l2)

           
println(L2)

# eval(VTK_mix_pressure)

dₚᵤ = zeros(2*nₚ + nᵤ)
dₚ = zeros(2*nₚ)
dᵤ = zeros(nᵤ)
for (i,node) in enumerate(nodes_p)
    x = node.x
    y = node.y
    dₚᵤ[2*i-1] = - ∂T∂x(x,y)
    dₚᵤ[2*i]   = - ∂T∂y(x,y)
    dₚ[2*i-1] = - ∂T∂x(x,y)
    dₚ[2*i]   = - ∂T∂y(x,y)
end
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    dₚᵤ[2*nₚ+i] = T(x,y)
    dᵤ[i] = T(x,y)
end

err1 = kₚₚ*dₚ - kₚᵤ'*dᵤ
err2 = kₚₙ*dᵤ - fₚ
err3 = -(kₚᵤ+kₚₙ')*dₚ-fᵤ
err4 = kₚₚ*dₚ -(kₚᵤ'+kₚₙ)*dᵤ + fₚ
err5 = k*dₚᵤ-f
err6 = k*d-f
