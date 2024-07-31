
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_plate_with_hole.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 4
# ndiv2 = 20
nₚ = 80
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes, nodes_u= import_patchtest_mix("./msh/plate_with_hole_new_"*string(ndiv)*".msh","./msh/plate_with_hole_new_"*string(ndiv2)*".msh")
elements,nodes, nodes_u= import_patchtest_mix("./msh/plate_with_hole_new_"*string(ndiv)*".msh","./msh/plate_with_hole_new_bubble_"*string(nₚ)*".msh")
nₚ = length(nodes)
nᵤ = length(nodes_u)
# nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set𝝭!(elements["Ωᵘ"])
set𝝭!(elements["∂Ωᵘ"])
set𝝭!(elements["Γ¹ᵗᵘ"])
set𝝭!(elements["Γ²ᵗᵘ"])
set𝝭!(elements["Γ¹ᵍᵘ"])
set𝝭!(elements["Γ²ᵍᵘ"])
set𝝭!(elements["Γ³ᵍᵘ"])
set∇𝝭!(elements["Ωᵖ"])
set𝝭!(elements["∂Ωᵖ"])
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
n =3
# T(x,y) = 0.5*a₀+a₁*(r(x,y)+R^2/r(x,y))*cos(θ(x,y))+a₂*(r(x,y)^2+R^4/r(x,y)^2)*cos(2*θ(x,y))+a₃*(r(x,y)^3+R^6/r(x,y)^3)*cos(3*θ(x,y))+a₄*(r(x,y)^4+R^8/r(x,y)^4)*cos(4*θ(x,y))
T(x,y) = 0.5*a₀+a₁*(r(x,y)+R^2/r(x,y))*cos(θ(x,y))
∂T∂x(x,y) = a₁*(x/(x^2+y^2)^0.5-x/(x^2+y^2)^1.5)/(1+y^2/x^2)^0.5+a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)*y^2/((1+y^2/x^2)^1.5*x^3)
∂T∂y(x,y) = a₁*(y/(x^2+y^2)^0.5-y/(x^2+y^2)^1.5)/(1+y^2/x^2)^0.5-a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)*y/((1+y^2/x^2)^1.5*x^2)
P₁(x,y) = -∂T∂x(x,y)
P₂(x,y) = -∂T∂y(x,y)
∂²T∂x²(x,y)  = a₁*(-x^2/(x^2+y^2)^1.5+1/(x^2+y^2)^0.5+3*x^2/(x^2+y^2)^2.5-1/(x^2+y^2)^1.5)/(1+y^2/x^2)^0.5 + 2*a₁*(x/(x^2+y^2)^0.5-x/(x^2+y^2)^1.5)*y^2/(1+y^2/x^2)^1.5/x^3+3*a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)*y^4/(1+y^2/x^2)^2.5/x^6-3*a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)*y^2/(1+y^2/x^2)^1.5/x^4
∂²T∂y²(x,y)  = a₁*(-y^2/(x^2+y^2)^1.5+1/(x^2+y^2)^0.5+3*y^2/(x^2+y^2)^2.5-1/(x^2+y^2)^1.5)/(1+y^2/x^2)^0.5 - 2*a₁*(y/(x^2+y^2)^0.5-y/(x^2+y^2)^1.5)*y/(1+y^2/x^2)^1.5/x^2+3*a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)*y^2/(1+y^2/x^2)^2.5/x^4-a₁*((x^2+y^2)^0.5+1/(x^2+y^2)^0.5)/(1+y^2/x^2)^1.5/x^2
b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))
# b(x,y) = 0.0 

# T(x,y) = (x+y)^n
# ∂T∂x(x,y) = n*(x+y)^abs(n-1)
# ∂T∂y(x,y) = n*(x+y)^abs(n-1)
# P₁(x,y) = -∂T∂x(x,y)
# P₂(x,y) = -∂T∂y(x,y)
# ∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
# b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))


eval(prescribe)
prescribe!(elements["∂Ωᵘ"],:g=>(x,y,z)->0.0)
ops = [
       Operator{:∫vtdΓ}(),
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
opsᵛ[3](elements["Γ¹ᵍᵖ"],elements["Γ¹ᵍᵘ"],kₚₙ,fₚ)
opsᵛ[3](elements["Γ²ᵍᵖ"],elements["Γ²ᵍᵘ"],kₚₙ,fₚ)
opsᵛ[3](elements["Γ³ᵍᵖ"],elements["Γ³ᵍᵘ"],kₚₙ,fₚ)
# ops[1](elements["Γ²ᵗᵘ"],fᵤ)
ops[1](elements["Γ¹ᵗᵘ"],fᵤ)
# ops[1](elements["Γ²ᵍᵘ"],fᵤ) 
ops[3](elements["Ωᵘ"],fᵤ)

# kₚᵤ = kₚᵤ+kₚₙ

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
err2 = kₚₙ*dᵤ + fₚ
err3 = (kₚᵤ+kₚₙ)'*dₚ-fᵤ
err4 = kₚₚ*dₚ - (kₚᵤ+kₚₙ)*dᵤ + fₚ
err5 = k*dₚᵤ-f
err6 = k*d-f
