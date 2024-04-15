using Revise, ApproxOperator, LinearAlgebra

include("import_patchtest.jl")
ndiv = 3

elements, nodes, fig = import_patchtest_tri6("./msh/patchtest_tri6_"*string(ndiv)*".msh")
# elements, nodes, fig = import_patchtest_stripe("./msh/patchtest_"*string(ndiv)*".msh")
# elements, nodes, fig = import_patchtest_unionJack("./msh/patchtest_"*string(ndiv)*".msh")
# elements, nodes, fig = import_patchtest_cross("./msh/patchtest_"*string(ndiv)*".msh")

nₚ = length(nodes)

Ē = 1.0
ν̄ = 0.3
# ν̄ = 0.49999
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

n = 10
u(x,y) = (x+y)^n
v(x,y) = -(x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂v∂x(x,y) = -n*(x+y)^abs(n-1)
∂v∂y(x,y) = -n*(x+y)^abs(n-1)
∂²u∂x²(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²v∂x²(x,y) = -n*(n-1)*(x+y)^abs(n-2)
∂²v∂x∂y(x,y) = -n*(n-1)*(x+y)^abs(n-2)
∂²v∂y²(x,y) = -n*(n-1)*(x+y)^abs(n-2)
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

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])
# set∇𝝭!(elements["Ωᵛ"])
# set𝝭!(elements["Γᵍ"])

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e10*E),
       Operator{:∫∫vᵢbᵢdxdy}(),
       Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
]
opsᵛ = [
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]

kᵛ = zeros(2*nₚ,2*nₚ)
kᵈ = zeros(2*nₚ,2*nₚ)
kᵍ = zeros(2*nₚ,2*nₚ) 
f = zeros(2*nₚ)


opsᵛ[1](elements["Ω"],kᵛ)
opsᵈ[1](elements["Ω"],kᵈ)
ops[3](elements["Γ"],kᵍ,f)
ops[4](elements["Ω"],f)

# ops[3](elements["Γᵍ"],kᵍ,f)

# vᵈ = eigvals(kᵈ+kᵍ)
# vᵛ = eigvals(kᵛ)
# v = eigvals(kᵛ,kᵈ+kᵍ)

# fig
k = kᵛ+kᵈ+kᵍ
d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ωᵍ"])
error = ops[5](elements["Ωᵍ"])