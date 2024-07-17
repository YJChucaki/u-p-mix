
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_patchtest.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv = 11
nₚ = 20
i = 20
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_bubble_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p,Ω = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_quad_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p ,Ω = import_patchtest_mix("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_quad8_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes = import_patchtest_Q4P1("./msh/patchtest_quad_"*string(ndiv)*".msh")
# elements,nodes = import_patchtest_Q4R1("./msh/patchtest_quad_"*string(ndiv)*".msh")
# elements,nodes,nodes_p = import_patchtest_T6P3("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
# elements,nodes = import_patchtest_Q8P3("./msh/patchtest_quad8_"*string(ndiv)*".msh")
# elements,nodes = import_patchtest_fem("./msh/patchtest_quad8_"*string(ndiv)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)
## for Q4P1 or Q4R1
# nₚ = length(elements["Ωᵖ"])
# for Q8P3
# nₚ = 3*length(elements["Ωᵖ"])

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γ"])
Ē = 1.0
ν̄ = 0.499999
# ν̄ = 0.3
# ν = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

n = 1
u(x,y) = (x+y)^n
v(x,y) = -(x+y)^n
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

# u(x,y) = (1+2*x+3*y)^n
# v(x,y) = (4+5*x+6*y)^n
# ∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
# ∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
# ∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
# ∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
# ∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

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

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν̄),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e12*E),
       Operator{:∫∫vᵢbᵢdxdy}(),
       Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄),
       Operator{:Hₑ_Incompressible}(:E=>Ē,:ν=>ν̄)
]
opsᵛ = [
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]

kᵅ = zeros(2*nᵤ,2*nᵤ)
fᵅ = zeros(2*nᵤ)
kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
kₚₚ = zeros(nₚ,nₚ)
fᵤ = zeros(2*nᵤ)


opsᵈ[1](elements["Ω"],kᵤᵤ)


opsᵛ[1](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
opsᵛ[2](elements["Ωᵖ"],kₚₚ)
ops[3](elements["Γ"],kᵅ,fᵅ)
ops[4](elements["Ω"],fᵤ)


# kᵈ = kᵤᵤ
# kᵛ = kᵤₚ*(kₚₚ\kᵤₚ')
k = [kᵤᵤ+kᵅ kₚᵤ';kₚᵤ kₚₚ]
f = [fᵤ+fᵅ;zeros(nₚ)]
# d = (kᵛ+kᵈ)\f
kᵈ = kᵤᵤ
kᵛ = -kₚᵤ'*(kₚₚ\kₚᵤ)
vᵈ = eigvals(kᵈ)
vᵛ = eigvals(kᵛ)
γ = eigvals(kᵛ,kᵈ)
# println(γ[2*nᵤ-nₚ+1])
d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
p  = d[2*nᵤ+1:end]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>p)

set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵍᵖ"])
# h1,l2= ops[6](elements["Ωᵍ"])
h1,l2,h1_dil,h1_dev,p_error= ops[5](elements["Ωᵍ"],elements["Ωᵍᵖ"])
L2 = log10(l2)
H1 = log10(h1)
H1_dil = log10(h1_dil)
H1_dev = log10(h1_dev)
P_error = log10(p_error)        
println(L2,H1)
println(H1_dil,H1_dev)
println(P_error)
# eval(VTK_mix_pressure)
# println(l2,h1)
# println(log10(sqrt(γ[1])))
# println(h1_dil,h1_dev)
# @save compress=true "jld/patchtest_mix_tri3_bubble_"*string(nₚ)*".jld" q


# d̄ = zeros(2*nᵤ+nₚ)
# d̃ = zeros(2*nᵤ)
# d̄₁ = zeros(nᵤ)
# d̄₂ = zeros(nᵤ)
# p̄ = zeros(nₚ)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     d̄₁[i] = u(x,y)
#     d̄₂[i] = v(x,y)
#     d̄[2*i-1] = u(x,y)
#     d̄[2*i] = v(x,y)
#     d̃[2*i-1] = u(x,y)
#     d̃[2*i] = v(x,y)
# end

# for (i,node) in enumerate(nodes_p)
#     x = node.x
#     y = node.y
#     p̄[i] = 𝑝(x,y)
#     d̄[2*nᵤ+i] = 𝑝(x,y)
# end

# err_d₁ = d₁ - d̄₁
# err_d₂ = d₂ - d̄₂
# err_p = p - p̄

# # err = k*d̄ .- f
# # err = kᵅ*d̃ .- fᵅ
# # err = [kᵤᵤ kᵤₚ]*d̄ - fᵤ
# # err = kᵤᵤ*d̃
# err = kᵤₚ*p̄