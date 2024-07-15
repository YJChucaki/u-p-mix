
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 11
nᵤ = 15

elements, nodes, nodes_u = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nᵤ)*".msh")

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
∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))

eval(prescribe)
prescribe!(elements["∂Ωᵘ"],:g=>(x,y,z)->0.0)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
       Operator{:∫vbdΩ}(),
       Operator{:L₂}(),
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
opsᵛ[3](elements["Γᵖ"],elements["Γᵘ"],kₚₙ,fₚ)
ops[3](elements["Ωᵘ"],fᵤ)

k = [kₚₚ -kₚᵤ-kₚₙ;-kₚᵤ'-kₚₙ' kᵤᵤ]
f = [-fₚ;-fᵤ]
d = k\f
p₁ = d[1:2:2*nₚ] 
p₂ = d[2:2:2*nₚ]
u  = d[2*nₚ+1:end]

push!(nodes,:d₁=>p₁,:d₂=>p₂)
push!(nodes_u,:d=>u)


set𝝭!(elements["Ωᵍᵘ"])
l2= ops[4](elements["Ωᵍᵘ"])
L2 = log10(l2)

           
println(L2)

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
# err5 = k*dₚᵤ-f
# err6 = k*d-f
