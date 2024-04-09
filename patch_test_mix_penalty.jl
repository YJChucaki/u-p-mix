
using ApproxOperator, Tensors, JLD

ndiv= 11
 i=105

include("import_prescrible_ops.jl")
include("import_patchtest.jl")
nₚ = 50
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)
 
Ē = 3e6
ν̄ = 0.4999999999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)
    
n = 3
u(x,y) = (x+y)^n
v(x,y) = -(x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂v∂x(x,y) = -n*(x+y)^abs(n-1)
∂v∂y(x,y) = -n*(x+y)^abs(n-1)

eval(prescribe)


set∇𝝭!(elements["Ωᵘ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Ωᵍᵖ"])
set𝝭!(elements["Γ¹"])
set𝝭!(elements["Γ²"])
set𝝭!(elements["Γ³"])
set𝝭!(elements["Γ⁴"])

eval(opsupmix)
kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kᵤₚ = zeros(2*nᵤ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(2*nᵤ)
fp= zeros(nₚ)
opsup[3](elements["Ωᵘ"],kᵤᵤ)
opsup[4](elements["Ωᵘ"],elements["Ωᵖ"],kᵤₚ)
opsup[5](elements["Ωᵖ"],kₚₚ)
# opsup[6](elements["Γᵗ"],f)

αᵥ = 1e9
eval(opsPenalty)

opsα[1](elements["Γ¹"],kᵤᵤ,f)
opsα[1](elements["Γ²"],kᵤᵤ,f)
opsα[1](elements["Γ³"],kᵤᵤ,f)
opsα[1](elements["Γ⁴"],kᵤᵤ,f)

# kᵈ = kᵤᵤ
# kᵛ = kᵤₚ*(kₚₚ\kᵤₚ')
k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;fp]
# d = (kᵛ+kᵈ)\f

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
q  = d[2*nᵤ+1:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>q)

# @save compress=true "jld/patchtest_gauss_penalty.jld" d₁ d₂ d₃

set∇𝝭!(elements["Ωᵍ"])
h1,l2,h1_dil,h1_dev= opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
# h1,l2,h1_dil,h1_dev = opsup[8](elements["Ω"],elements["Ωᵖ"])
L2 = log10(l2)
H1 = log10(h1)
H1_dil = log10(h1_dil)
H1_dev = log10(h1_dev)
# h = log10(12.0/ndiv)
           
# println(L2,H1)
println(l2,h1)
# println(H1_dil,H1_dev)