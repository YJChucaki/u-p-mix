using Revise, ApproxOperator, LinearAlgebra

include("import_patchtest.jl")
ndiv= 3
nₚ = 100
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")

nᵤ = length(nodes)
nₚ = length(nodes_p)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γ"])

Ē = 1.0
ν̄ = 0.3
# ν̄ = 0.499999
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

n = 10
u(x,y) = (x+y)^n
v(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂v∂x(x,y) = n*(x+y)^abs(n-1)
∂v∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²v∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²v∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²v∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
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

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e10*E),
       Operator{:∫∫vᵢbᵢdxdy}(),
       Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄)
]
opsᵛ = [
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kᵤₚ = zeros(2*nᵤ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(2*nᵤ)

# opsᵈ[1](elements["Ω"],kᵤᵤ)
# opsᵛ[1](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
# opsᵛ[2](elements["Ωᵖ"],kₚₚ)
# ops[3](elements["Γ"],kᵤᵤ,f)
ops[4](elements["Ω"],f)

# kᵈ = kᵤᵤ
# kᵛ = -kᵤₚ*(kₚₚ\kᵤₚ')
# vᵈ = eigvals(kᵈ)
# vᵛ = eigvals(kᵛ)
# v = eigvals(kᵛ,kᵈ)

# d = (kᵛ+kᵈ)\f
# fig

# k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
# f = [f;zeros(nₚ)]

k = zeros(2*nᵤ,2*nᵤ)
ops[1](elements["Ω"],k)
ops[3](elements["Γ"],k,f)

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
# p = d[2*nᵤ+1:end]
p = zeros(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p, :q=>p)

set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵍᵖ"])
# error = ops[5](elements["Ωᵍ"], elements["Ωᵍᵖ"])
h1,l2,h1_dil,h1_dev= ops[5](elements["Ωᵍ"],elements["Ωᵍᵖ"])
L2 = log10(l2)
H1 = log10(h1)
H1_dil = log10(h1_dil)
H1_dev = log10(h1_dev)
           
# println(L2,H1)
println(l2,h1)
# println(H1_dil,H1_dev)
