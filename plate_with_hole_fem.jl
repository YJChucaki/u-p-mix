
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_plate_with_hole.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 16
# ndiv2 =4
# nₚ =140
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes= import_patchtest_fem("./msh/plate_with_hole_quad_"*string(ndiv)*".msh")
elements,nodes= import_patchtest_Q4P1("./msh/plate_with_hole_quad_"*string(ndiv)*".msh")
# elements,nodes= import_patchtest_fem("./msh/plate_with_hole_new_quad_"*string(ndiv)*".msh")
# elements,nodes= import_patchtest_fem("./msh/plate_with_hole.msh")
nᵤ = length(nodes)
# nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set∇𝝭!(elements["Ω"])

set𝝭!(elements["Γ¹ᵗ"])
set𝝭!(elements["Γ²ᵗ"])
set𝝭!(elements["Γ¹ᵍ"])
set𝝭!(elements["Γ²ᵍ"])
set𝝭!(elements["Γ³ᵍ"])
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
# n = 1
# T(x,y) = (x+y)^n
# ∂T∂x(x,y) = n*(x+y)^abs(n-1)
# ∂T∂y(x,y) = n*(x+y)^abs(n-1)
# P₁(x,y) = -∂T∂x(x,y)
# P₂(x,y) = -∂T∂y(x,y)
# ∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
# b(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))


eval(prescribeForFem)

ops = [
       Operator{:∫vtdΓ}(),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
       Operator{:∫vbdΩ}(),
       Operator{:L₂}(),
       Operator{:L₂_heat_flux}(),
       Operator{:∫∫∇TᵢD∇Tⱼdxdy}(:D=>D,:t=>t),
]


kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
k = zeros(nᵤ,nᵤ)
f = zeros(nᵤ)


ops[6](elements["Ω"],k)
ops[3](elements["Ω"],f)
ops[2](elements["Γ¹ᵍ"],kᵅ,f)
ops[2](elements["Γ²ᵍ"],kᵅ,f)
ops[2](elements["Γ³ᵍ"],kᵅ,f)
ops[2](elements["Γ¹ᵗ"],kᵅ,f)
# ops[2](elements["Γ²ᵗ"],kᵅ,f)
# ops[1](elements["Γ²ᵍ"],f)
d = (k+kᵅ)\f #temperatures
push!(nodes,:d=>d)
# p₁ = zeros(nᵤ)
# p₂ = zeros(nᵤ)
# for ap in elements["Ω"]
#        𝓒 = ap.𝓒
#        𝓖 = ap.𝓖
       
#        for (i,ξ) in enumerate(𝓖)
#                B₁ = ξ[:∂𝝭∂x]
#                B₂ = ξ[:∂𝝭∂y]
#                for (j,xⱼ) in enumerate(𝓒)
#                    I = xⱼ.𝐼
#                    p₁[I] -= B₁[j]*xⱼ.d
#                    p₂[I] -= B₂[j]*xⱼ.d
#                end 
#        end
#    end
#    push!(nodes,:d₁=>p₁,:d₂=>p₂)


set∇𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍᵖ"])
l2_u= ops[4](elements["Ωᵍ"])
l2_p= ops[5](elements["Ωᵍ"])

L2_u = log10(l2_u)
L2_p = log10(l2_p)
     
println(L2_u)
println(L2_p)
           


# eval(VTK_mix_pressure)

# dₚᵤ = zeros(2*nₚ + nᵤ)
# dₚ = zeros(2*nₚ)
# dᵤ = zeros(nᵤ)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     dₚᵤ[2*i-1] = -∂T∂x(x,y)
#     dₚᵤ[2*i]   = -∂T∂y(x,y)
#     dₚ[2*i-1] = -∂T∂x(x,y)
#     dₚ[2*i]   = -∂T∂y(x,y)
# end
# for (i,node) in enumerate(nodes_u)
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
