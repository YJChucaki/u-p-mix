using Revise, ApproxOperator, LinearAlgebra

include("import_patchtest.jl")
ndiv = 21

# elements, nodes, fig = import_patchtest_stripe("./msh/cantilever_8.msh")
# elements, nodes, fig = import_patchtest_stripe("./msh/patchtest_"*string(ndiv)*".msh")
# elements, nodes, fig = import_patchtest_unionJack("./msh/patchtest_"*string(ndiv)*".msh")
elements, nodes, fig = import_patchtest_cross("./msh/patchtest_"*string(ndiv)*".msh")

nₚ = length(nodes)

Ē = 1.0
ν̄ = 0.49999
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Ωᵛ"])
# set∇𝝭!(elements["Ωᵛ"])
# set𝝭!(elements["Γᵍ"])

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e9*E),
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
# ops[3](elements["Γᵍ"],kᵍ,f)

# ops[3](elements["Γᵍ"],kᵍ,f)

vᵈ = eigvals(kᵈ+kᵍ)
vᵛ = eigvals(kᵛ)
v = eigvals(kᵛ,kᵈ+kᵍ)

fig