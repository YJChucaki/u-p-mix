
using ApproxOperator, LinearAlgebra

include("import_patchtest_bubble.jl")

ndiv = 11
nₚ = 200
elements, nodes, nodes_p = import_mix_bubble("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
nᵤ = length(nodes)
nₒ = length(elements["Ω"])
nₛ = length(elements["Ω"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωᵇ"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γᵖ"])
set𝝭!(elements["Γˢ"])

Ē = 1.0
ν̄ = 0.4999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

n = 3
u(x,y) = (1+2*x+3*y)^n
v(x,y) = (4+5*x+6*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = Ē/(1+ν̄)/(1-2*ν̄)*((1-ν̄)*ε₁₁(x,y) + ν̄*ε₂₂(x,y))
σ₂₂(x,y) = Ē/(1+ν̄)/(1-2*ν̄)*(ν̄*ε₁₁(x,y) + (1-ν̄)*ε₂₂(x,y))
σ₃₃(x,y) = Ē*ν̄/(1+ν̄)/(1-2*ν̄)*(ε₁₁(x,y) + ε₂₂(x,y))
σ₁₂(x,y) = Ē/(1+ν̄)*ε₁₂(x,y)
𝑝(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3
𝑠₁₁(x,y) = Ē/(1+ν̄)*( 2/3*ε₁₁(x,y) - 1/3*ε₂₂(x,y))
𝑠₂₂(x,y) = Ē/(1+ν̄)*(-1/3*ε₁₁(x,y) + 2/3*ε₂₂(x,y))
𝑠₁₂(x,y) = Ē/(1+ν̄)*ε₁₂(x,y)
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))
∂σ₁₁∂x(x,y) = E/(1-ν^2)*(∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1-ν^2)*(∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂x(x,y) + ∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂y(x,y) + ∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)

eval(prescribe)

opsᵖ = [
    Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫pnᵢgᵢds}(),
]

opsˢ = [
    Operator{:∫∫δsᵢⱼsᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫∫sᵢⱼεᵢⱼdxdy}(),
    Operator{:∫sᵢⱼnⱼgᵢds}(),
]

ops = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢgᵢds}(:α=>1e15*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
]

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
kₛᵤ = zeros(4*nₛ,2*nᵤ)
kₒᵤ = zeros(2*nₒ,2*nᵤ)
kₚₚ = zeros(nₚ,nₚ)
kₛₛ = zeros(4*nₛ,4*nₛ)
kₒₒ = zeros(2*nₒ,2*nₒ)
fᵤ = zeros(2*nᵤ)
fₚ = zeros(nₚ)
fₛ = zeros(4*nₛ)
fₒ = zeros(2*nₒ)

opsᵖ[1](elements["Ωᵖ"],kₚₚ)
opsᵖ[2](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
# opsᵖ[3](elements["Γ"],elements["Γᵖ"],kₚᵤ,fₚ)

# opsˢ[1](elements["Ωˢ"],kₛₛ)
# opsˢ[2](elements["Ω"],elements["Ωˢ"],kₛᵤ)
# opsˢ[3](elements["Γ"],elements["Γˢ"],kₛᵤ,fₛ)

ops[1](elements["Ω"],kᵤᵤ)
ops[1](elements["Ωᵇ"],kₒₒ)
ops[3](elements["Ω"],fᵤ)
ops[3](elements["Ωᵇ"],fₒ)
ops[4](elements["Γ"],kᵤᵤ,fᵤ)


# k = [kᵤᵤ kₚᵤ' kₛᵤ' kₒᵤ';
#      kₚᵤ kₚₚ zeros(nₚ,4*nₛ) zeros(nₚ,2*nₒ);
#      kₛᵤ zeros(4*nₛ,nₚ) kₛₛ zeros(4*nₛ,2*nₒ);
#      kₒᵤ zeros(2*nₒ,nₚ) zeros(2*nₒ,4*nₛ) kₒₒ]
k = [kᵤᵤ kₚᵤ' kₒᵤ';
     kₚᵤ kₚₚ zeros(nₚ,2*nₒ);
     kₒᵤ zeros(2*nₒ,nₚ) kₒₒ]
f = [fᵤ;fₚ;fₒ]

# # k = zeros(2*nᵤ,2*nᵤ)
# # f = zeros(2*nᵤ)
# # ops[2](elements["Ω"],f)
# # ops[4](elements["Ω"],k)
# # ops[5](elements["Γ"],k,f)

d = k\f

# dᵤ = zeros(2*nᵤ)
# dₚ = zeros(nₚ)
# dₛ = zeros(4*nₛ)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     dᵤ[2*i-1] = u(x,y)
#     dᵤ[2*i]   = v(x,y)
# end

# for (i,node) in enumerate(nodes_p)
#     x = node.x
#     y = node.y
#     dₚ[i] = 𝑝(x,y)
# end

# for (i,element) in enumerate(elements["Ω"])

#     # dₛ[12*i-11] = σ₁₁(0,0)
#     # dₛ[12*i-10] = σ₂₂(0,0)
#     # dₛ[12*i-9]  = σ₃₃(0,0)
#     # dₛ[12*i-8]  = σ₁₂(0,0)
#     # dₛ[12*i-7]  = σ₁₁(1,0) - σ₁₁(0,0)
#     # dₛ[12*i-6]  = σ₂₂(1,0) - σ₂₂(0,0)
#     # dₛ[12*i-5]  = σ₃₃(1,0) - σ₃₃(0,0)
#     # dₛ[12*i-4]  = σ₁₂(1,0) - σ₁₂(0,0)
#     # dₛ[12*i-3]  = σ₁₁(0,1) - σ₁₁(0,0)
#     # dₛ[12*i-2]  = σ₂₂(0,1) - σ₂₂(0,0)
#     # dₛ[12*i-1]  = σ₃₃(0,1) - σ₃₃(0,0)
#     # dₛ[12*i]    = σ₁₂(0,1) - σ₁₂(0,0)

#     dₛ[4*i-3] = σ₁₁(0,0)
#     dₛ[4*i-2] = σ₂₂(0,0)
#     dₛ[4*i-1] = σ₃₃(0,0)
#     dₛ[4*i] = σ₁₂(0,0)
#     dₛ[3*i-2] = E/(1+ν)*( 2/3*2 - 1/3*6)
#     dₛ[3*i-1] = E/(1+ν)*(-1/3*2 + 2/3*6)
#     dₛ[3*i]   = E/(1+ν)*0.5*(3+5)
# end
# # dₛ .= d[2*nᵤ+nₚ+1:end]

# # errₚ = kₚₚ*dₚ + kₚᵤ*dᵤ - fₚ
# # errₛ = kₛₛ*dₛ + kₛᵤ*dᵤ - fₛ
# # err = kₛₛ*dₛ + kₛᵤ*dᵤ
# # errᵤ = kₚᵤ'*dₚ + kₛᵤ'*dₛ - fᵤ
# # err = kₚᵤ'*dₚ
# # err = kₛᵤ'*dₛ
# # norm(errᵤ)
# # norm(errₚ)
# # norm(errₛ)

d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ωᵍ"])
h1,l2= ops[5](elements["Ωᵍ"])
L2 = log10(l2)
H1 = log10(h1)
           
println(L2)
println(H1)