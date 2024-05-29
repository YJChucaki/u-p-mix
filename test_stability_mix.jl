
using ApproxOperator, Tensors,  LinearAlgebra, JLD
include("import_patchtest.jl")
# for i=2:10
   
ndiv= 9
nₚ = 90
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_bubble_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_quad_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_quad8_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
# elements,nodes = import_patchtest_Q4P1("./msh/patchtest_quad_"*string(ndiv)*".msh")
# elements,nodes = import_patchtest_Q4R1("./msh/patchtest_quad_"*string(ndiv)*".msh")
# elements,nodes,nodes_p = import_patchtest_T6P3("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
# elements,nodes = import_patchtest_Q8P3("./msh/patchtest_quad8_"*string(ndiv)*".msh")
nᵤ = length(nodes)
# nₚ = length(nodes_p)

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γ"])
Ē = 1.0
ν̄ = 0.4999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:t₁=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:t₂=>(x,y,z)->0.0)

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e13*E),
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

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kᵤₚ = zeros(2*nᵤ,nₚ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(2*nᵤ)


opsᵈ[1](elements["Ω"],kᵤᵤ)


opsᵛ[1](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
opsᵛ[2](elements["Ωᵖ"],kₚₚ)
ops[3](elements["Γ¹"],kᵤᵤ,f)
ops[3](elements["Γ⁴"],kᵤᵤ,f)
ops[2](elements["Γ²"],f)


k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;zeros(nₚ)]
d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
p  = d[2*nᵤ+1:end]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>p)

@save compress=true "jld/stability_tri3_"*string(nₚ)*".jld" d₁ d₂ p
