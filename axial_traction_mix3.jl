
using ApproxOperator, LinearAlgebra

include("import_patchtest_mix3.jl")

ndiv = 9
nₚ = 140
elements, nodes, nodes_p = import_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
nᵤ = length(nodes)
nₛ = length(elements["Ω"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γᵖ"])
set𝝭!(elements["Γˢ"])

P = 1.0
Ē = 1.0
ν̄ = 0.499999999999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

u(x,y) = (1-ν̄^2)/Ē*P*x
v(x,y) = -ν̄*(1+ν̄)/Ē*P*y
∂u∂x(x,y) = (1-ν̄^2)/Ē*P
∂u∂y(x,y) = 0.0
∂v∂x(x,y) = 0.0
∂v∂y(x,y) = -ν̄*(1+ν̄)/Ē*P
∂²u∂x²(x,y)  = 0.0
∂²u∂x∂y(x,y) = 0.0
∂²u∂y²(x,y)  = 0.0
∂²v∂x²(x,y)  = 0.0
∂²v∂x∂y(x,y) = 0.0
∂²v∂y²(x,y)  = 0.0

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

eval(prescribe_at)

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
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e10*E)
]

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
kₛᵤ = zeros(4*nₛ,2*nᵤ)
kₚₚ = zeros(nₚ,nₚ)
kₛₛ = zeros(4*nₛ,4*nₛ)
fᵤ = zeros(2*nᵤ)
fₚ = zeros(nₚ)
fₛ = zeros(4*nₛ)

opsᵖ[1](elements["Ωᵖ"],kₚₚ)
opsᵖ[2](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
opsᵖ[3](elements["Γ¹"],elements["Γ¹ᵖ"],kₚᵤ,fₚ)
opsᵖ[3](elements["Γ⁴"],elements["Γ⁴ᵖ"],kₚᵤ,fₚ)

opsˢ[1](elements["Ωˢ"],kₛₛ)
opsˢ[2](elements["Ω"],elements["Ωˢ"],kₛᵤ)
opsˢ[3](elements["Γ¹"],elements["Γ¹ˢ"],kₛᵤ,fₛ)
opsˢ[3](elements["Γ⁴"],elements["Γ⁴ˢ"],kₛᵤ,fₛ)

ops[1](elements["Γ²"],fᵤ)

k = [zeros(2*nᵤ,2*nᵤ) kₚᵤ' kₛᵤ';
     kₚᵤ kₚₚ zeros(nₚ,4*nₛ);
     kₛᵤ zeros(4*nₛ,nₚ) kₛₛ]
f = [fᵤ;fₚ;fₛ]

# k = zeros(2*nᵤ,2*nᵤ)
# f = zeros(2*nᵤ)
# ops[2](elements["Ω"],f)
# ops[4](elements["Ω"],k)
# ops[5](elements["Γ"],k,f)

d = k\f

dᵤ = zeros(2*nᵤ)
dₚ = zeros(nₚ)
dₛ = zeros(4*nₛ)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    dᵤ[2*i-1] = u(x,y)
    dᵤ[2*i]   = v(x,y)
end

for (i,node) in enumerate(nodes_p)
    x = node.x
    y = node.y
    dₚ[i] = 𝑝(x,y)
end

for (i,element) in enumerate(elements["Ω"])

    # dₛ[12*i-11] = σ₁₁(0,0)
    # dₛ[12*i-10] = σ₂₂(0,0)
    # dₛ[12*i-9]  = σ₃₃(0,0)
    # dₛ[12*i-8]  = σ₁₂(0,0)
    # dₛ[12*i-7]  = σ₁₁(1,0) - σ₁₁(0,0)
    # dₛ[12*i-6]  = σ₂₂(1,0) - σ₂₂(0,0)
    # dₛ[12*i-5]  = σ₃₃(1,0) - σ₃₃(0,0)
    # dₛ[12*i-4]  = σ₁₂(1,0) - σ₁₂(0,0)
    # dₛ[12*i-3]  = σ₁₁(0,1) - σ₁₁(0,0)
    # dₛ[12*i-2]  = σ₂₂(0,1) - σ₂₂(0,0)
    # dₛ[12*i-1]  = σ₃₃(0,1) - σ₃₃(0,0)
    # dₛ[12*i]    = σ₁₂(0,1) - σ₁₂(0,0)

    dₛ[4*i-3] = σ₁₁(0,0)
    dₛ[4*i-2] = σ₂₂(0,0)
    dₛ[4*i-1] = σ₃₃(0,0)
    dₛ[4*i] = σ₁₂(0,0)
    dₛ[3*i-2] = E/(1+ν)*( 2/3*2 - 1/3*6)
    dₛ[3*i-1] = E/(1+ν)*(-1/3*2 + 2/3*6)
    dₛ[3*i]   = E/(1+ν)*0.5*(3+5)
end
# dₛ .= d[2*nᵤ+nₚ+1:end]

# errₚ = kₚₚ*dₚ + kₚᵤ*dᵤ - fₚ
# errₛ = kₛₛ*dₛ + kₛᵤ*dᵤ - fₛ
# err = kₛₛ*dₛ + kₛᵤ*dᵤ
# errᵤ = kₚᵤ'*dₚ + kₛᵤ'*dₛ - fᵤ
# err = kₚᵤ'*dₚ
# err = kₛᵤ'*dₛ
# norm(errᵤ)
# norm(errₚ)
# norm(errₛ)

d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ωᵍ"])
h1,l2= ops[3](elements["Ωᵍ"][1:1])
L2 = log10(l2)
H1 = log10(h1)
           
println(L2)
println(H1)