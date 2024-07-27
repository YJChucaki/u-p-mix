
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_plate_with_hole.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 6
ndiv2 =3
# nₚ = 15
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
elements,nodes,nodes_p= import_patchtest_mix_old("./msh/plate_with_hole_"*string(ndiv)*".msh","./msh/plate_with_hole_"*string(ndiv2)*".msh")

nₚ = length(nodes_p)
nᵤ = length(nodes)
# nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set∇𝝭!(elements["Ωᵘ"])

set𝝭!(elements["Γ¹ᵗᵘ"])
set𝝭!(elements["Γ²ᵗᵘ"])
set𝝭!(elements["Γ¹ᵍᵘ"])
set𝝭!(elements["Γ²ᵍᵘ"])
set𝝭!(elements["Γ³ᵍᵘ"])
set∇𝝭!(elements["Ωᵖ"])

set𝝭!(elements["Γ¹ᵗᵖ"])
set𝝭!(elements["Γ²ᵗᵖ"])
set𝝭!(elements["Γ¹ᵍᵖ"])
set𝝭!(elements["Γ²ᵍᵖ"])
set𝝭!(elements["Γ³ᵍᵖ"])
D=1   #thermal conductivity coefficient
t=1 #thickness
R = 1.0
r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
a₀ = 249.9977
a₁ = 24.23894
a₂ = 0.000025
a₃ = -0.141899
a₄ = 0.00007
n =1
# T(x,y) = 0.5*a₀+a₁*(r(x,y)+R^2/r(x,y))*cos(θ(x,y))+a₂*(r(x,y)^2+R^4/r(x,y)^2)*cos(2*θ(x,y))+a₃*(r(x,y)^3+R^6/r(x,y)^3)*cos(3*θ(x,y))+a₄*(r(x,y)^4+R^8/r(x,y)^4)*cos(4*θ(x,y))
# T(x,y) = 0.5*a₀+a₁*(r(x,y)+R^2/r(x,y))*cos(θ(x,y))
# ∂T∂x(x,y) = a₁*(2*x-2*x/(x^2+y^2)^2)/(1+y^2/x^2)^0.5+a₁*(x^2+y^2+1/(x^2+y^2))*y^2/((1+y^2/x^2)^1.5*x^3)
# ∂T∂y(x,y) = a₁*(2*y-2*y/(x^2+y^2)^2)/(1+y^2/x^2)^0.5+a₁*(x^2+y^2+1/(x^2+y^2))*y/((1+y^2/x^2)^1.5*x^2)
T(x,y) = (x+y)^n
∂T∂x(x,y) = n*(x+y)^abs(n-1)
∂T∂y(x,y) = n*(x+y)^abs(n-1)
P₁(x,y) = -∂T∂x(x,y)
P₂(x,y) = -∂T∂y(x,y)
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
opsᵛ[2](elements["Γ¹ᵍᵖ"],elements["Γ¹ᵍᵘ"],kₚₙ,fₚ)
opsᵛ[2](elements["Γ²ᵍᵖ"],elements["Γ²ᵍᵘ"],kₚₙ,fₚ)
opsᵛ[2](elements["Γ³ᵍᵖ"],elements["Γ³ᵍᵘ"],kₚₙ,fₚ)
opsᵛ[2](elements["Γ¹ᵗᵖ"],elements["Γ¹ᵗᵘ"],kₚₙ,fₚ)
opsᵛ[2](elements["Γ²ᵗᵖ"],elements["Γ²ᵗᵘ"],kₚₙ,fₚ)
ops[3](elements["Ωᵘ"],fᵤ)



# kₚᵤ = kₚᵤ+kₚₙ

k = [kₚₚ -kₚᵤ'-kₚₙ;-kₚᵤ-kₚₙ' kᵤᵤ]

# k = [kₚₚ -kₚᵤ;-kₚᵤ' kᵤᵤ]
f = [fₚ;-fᵤ]
d = k\f
p₁ = d[1:2:2*nₚ] 
p₂ = d[2:2:2*nₚ]
u  = d[2*nₚ+1:end]

push!(nodes_p,:d₁=>p₁,:d₂=>p₂)
push!(nodes,:d=>u)

set𝝭!(elements["Ωᵍᵘ"])
set𝝭!(elements["Ωᵍᵖ"])
l2_u= ops[4](elements["Ωᵍᵘ"])
# l2_p= ops[5](elements["Ωᵍᵖ"])
L2_u = log10(l2_u)
# L2_p = log10(l2_p)

           
println(L2_u)
# println(L2_p)
           


# eval(VTK_mix_pressure)

# dₚᵤ = zeros(2*nₚ + nᵤ)
# dₚ = zeros(2*nₚ)
# dᵤ = zeros(nᵤ)
# for (i,node) in enumerate(nodes_p)
#     x = node.x
#     y = node.y
#     dₚᵤ[2*i-1] = -∂T∂x(x,y)
#     dₚᵤ[2*i]   = -∂T∂y(x,y)
#     dₚ[2*i-1] = -∂T∂x(x,y)
#     dₚ[2*i]   = -∂T∂y(x,y)
# end
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     dₚᵤ[2*nₚ+i] = T(x,y)
#     dᵤ[i] = T(x,y)
# end

# err1 = kₚₚ*dₚ - kₚᵤ*dᵤ
# err2 = kₚₙ*dᵤ + fₚ
# err3 = (kₚᵤ+kₚₙ)'*dₚ-fᵤ
# err4 = kₚₚ*dₚ - (kₚᵤ+kₚₙ)*dᵤ + fₚ
# err5 = k*dₚᵤ-f
# err6 = k*d-f
