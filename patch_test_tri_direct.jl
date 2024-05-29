using Revise, ApproxOperator, LinearAlgebra, Printf
include("import_patchtest.jl")
ndiv=5
elements,nodes = import_patchtest_tri6("./msh/patchtest_"*string(ndiv)*".msh")


nₚ = length(nodes)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])

E = 3e6
ν=0.3
n = 2
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

# ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1-ν)*n₁+E/(1+ν)*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)*n₁+E/(1-ν)*n₂)

eval(prescribe)
ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:g₂}(),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
]

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
for ap in elements["Γ"]
    x, = ap.𝓒
    x.d₁ = u(x.x,x.y)
    x.d₂ = v(x.x,x.y)
end

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](elements["Ω"],k)
# ops[3](elements["Γᵗ"],f)
ops[4].(elements["Γ"],k=k,f=f,dof=:d₁)
ops[4].(elements["Γ"],k=k,f=f,dof=:d₂)

d = k\f
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ωᵍ"])
h1,l2 = ops[5](elements["Ωᵍ"])
L2 = log10(l2)
H1 = log10(h1)
println(L2,H1)