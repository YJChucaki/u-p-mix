
using ApproxOperator, Tensors, JLD

# ndiv= 11
 ndiv_p= 11
 i=80
# 40,60-3
# 80-4
# 100,120-5
# 160,200-7

include("import_prescrible_ops.jl")
include("import_patch_test.jl")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
# elements, nodes ,nodes_p = import_cantilever_mix_quad4("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix_tri3("./msh/patchtest.msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix_quad4("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
    nᵤ = length(nodes)
    nₚ = length(nodes_p)
    nₘ=21
    s =1.5*12/ndiv_p*ones(nₚ)
    P = 1000
    Ē = 3e6
    ν̄ = 0.4999999
    # ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 48
    D = 12
    I = D^3/12
    EI = E*I
    K=Ē/3/(1-2ν̄ )
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    u(x,y) = x+y
    ∂u∂x(x,y) = 1.0
    ∂u∂y(x,y) = 1.0
    v(x,y) = x+y
    ∂v∂x(x,y) = 1.0
    ∂v∂y(x,y) = 1.0

   eval(prescribeForGauss)
   eval(prescribeForPenalty)
#    eval(prescribeVariables)
#    eval(opsGauss)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
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
opsup[3](elements["Ω"],kᵤᵤ)
opsup[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
opsup[5](elements["Ωᵖ"],kₚₚ)
# opsup[6](elements["Γᵗ"],f)

αᵥ = 1e9
eval(opsPenalty)

opsα[1](elements["Γ¹"],kᵤᵤ,f)
opsα[1](elements["Γ²"],kᵤᵤ,f)
opsα[1](elements["Γ³"],kᵤᵤ,f)
opsα[1](elements["Γ⁴"],kᵤᵤ,f)


k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;fp]

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
q  = d[2*nᵤ+1:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>q)

# @save compress=true "jld/patchtest_gauss_penalty.jld" d₁ d₂ d₃

set∇𝝭!(elements["Ωᵍ"])
h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
# h1,l2,h1_dil,h1_dev = opsup[8](elements["Ω"],elements["Ωᵖ"])
L2 = log10(l2)
H1 = log10(h1)
H1_dil = log10(h1_dil)
H1_dev = log10(h1_dev)
# h = log10(12.0/ndiv)
println(L2,H1)
println(H1_dil,H1_dev)