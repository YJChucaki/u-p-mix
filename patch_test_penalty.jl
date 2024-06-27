
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_patchtest.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv= 11
# elements,nodes = import_patchtest_fem("./msh/patchtest_quad8_"*string(ndiv)*".msh")
# elements, nodes = import_patchtest_fem("./msh/patchtest_tri6_"*string(ndiv)*".msh")
# elements, nodes = import_patchtest_cross("./msh/patchtest_"*string(ndiv)*".msh")
# elements, nodes = import_patchtest_stripe("./msh/patchtest_"*string(ndiv)*".msh")
elements, nodes = import_patchtest_fem("./msh/patchtest_"*string(ndiv)*".msh")

nᵤ = length(nodes)
nₑ = length(elements["Ω"])


set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])
E = 1.0
ν= 0.4999999
# ν = 0.3
# E = Ē/(1.0-ν̄^2)
# ν = ν̄/(1.0-ν̄)

# n = 2
# u(x,y) = (x+y)^n
# v(x,y) = (x+y)^n
# ∂u∂x(x,y) = n*(x+y)^abs(n-1)
# ∂u∂y(x,y) = n*(x+y)^abs(n-1)
# ∂v∂x(x,y) = n*(x+y)^abs(n-1)
# ∂v∂y(x,y) = n*(x+y)^abs(n-1)
# ∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
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

eval(prescribe)

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν ),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e15*E),
       Operator{:∫∫vᵢbᵢdxdy}(),
       Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν)
]

kᵅ = zeros(2*nᵤ,2*nᵤ)
fᵅ = zeros(2*nᵤ)
k = zeros(2*nᵤ,2*nᵤ)
f = zeros(2*nᵤ)



ops[1](elements["Ω"],k)
ops[3](elements["Γ"],kᵅ,fᵅ)
ops[4](elements["Ω"],f)

f = f+fᵅ

d = (k+kᵅ)\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]


push!(nodes,:d₁=>d₁,:d₂=>d₂)




set∇𝝭!(elements["Ωᵍ"])
h1,l2= ops[5](elements["Ωᵍ"])
L2 = log10(l2)
H1 = log10(h1)
# H1_dil = log10(h1_dil)
# H1_dev = log10(h1_dev)
           
println(L2,H1)

# eval(VTK_mix_pressure)
# println(l2,h1)
# println(log10(sqrt(γ[1])))
# println(h1_dil,h1_dev)
# @save compress=true "jld/patchtest_mix_tri3_bubble_"*string(nₚ)*".jld" q
