using Revise, ApproxOperator, LinearAlgebra

include("import_patchtest.jl")
ndiv= 11
nₚ = 50
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")

nᵤ = length(nodes)
nₚ = length(nodes_p)

set∇𝝭!(elements["Ωᵘ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γ"])

Ē = 1.0
ν̄ = 0.4999999999999
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

n = 1
u(x,y) = (x+y)^n
v(x,y) = -(x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂v∂x(x,y) = -n*(x+y)^abs(n-1)
∂v∂y(x,y) = -n*(x+y)^abs(n-1)

eval(prescribe)

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e9*E),
       Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
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


opsᵈ[1](elements["Ωᵘ"],kᵤᵤ)
opsᵛ[1](elements["Ωᵘ"],elements["Ωᵖ"],kᵤₚ)
opsᵛ[2](elements["Ωᵖ"],kₚₚ)
ops[3](elements["Γ"],kᵤᵤ,f)

kᵈ = kᵤᵤ
kᵛ = kᵤₚ*(kₚₚ\kᵤₚ')
vᵈ = eigvals(kᵈ)
vᵛ = eigvals(kᵛ)
# v = eigvals(kᵛ,kᵈ)

# fig

# k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
# f = [f;zeros(nₚ)]

# d = k\f
d = (kᵛ+kᵈ)\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

set∇𝝭!(elements["Ωᵍ"])
Hₑ_PlaneStress = ops[4](elements["Ωᵍ"])
