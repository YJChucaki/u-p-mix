using ApproxOperator, Tensors, JLD

include("import_patch_test.jl")
include("import_prescrible_ops.jl")
# elements, nodes = import_patchtest_fem("./msh/patchtest.msh")
elements, nodes = import_patchtest_mf("./msh/patchtest.msh")

nᵤ = length(nodes)
nₘ=21
#mf
s = 2.5/10*ones(nᵤ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

E = 3e6
ν=0.3
# ν=0.49999999999999
Ē = E/(1-ν^2)
ν̄ = ν/(1-ν)

# u(x,y) =  2*x*y+x^2+y^2
# v(x,y) = -2*x*y-x^2-y^2
# ∂u∂x(x,y) = 2*x+2*y
# ∂u∂y(x,y) = 2*x+2*y
# ∂v∂x(x,y) = -2*x-2*y
# ∂v∂y(x,y) = -2*x-2*y
# ∂²u∂x²(x,y) = 2.0
# ∂²u∂x∂y(x,y) = 2.0
# ∂²u∂y²(x,y) = 2.0
# ∂²v∂x²(x,y) = -2.0
# ∂²v∂x∂y(x,y) = -2.0
# ∂²v∂y²(x,y) = -2.0

u(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
v(x,y) = x+y
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0
eval(prescribeForGauss)
eval(prescribeForPenalty)



set∇𝝭!(elements["Ωᵘ"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γ¹"])
set𝝭!(elements["Γ²"])
set𝝭!(elements["Γ³"])
set𝝭!(elements["Γ⁴"])

eval(opsFEM)
k = zeros(2*nᵤ,2*nᵤ)
kα = zeros(2*nᵤ,2*nᵤ)
f = zeros(2*nᵤ)
opsFEM[1](elements["Ω"],k)

αᵥ = 1e9
eval(opsPenalty)

opsα[1](elements["Γ¹"],kα,f)
opsα[1](elements["Γ²"],kα,f)
opsα[1](elements["Γ³"],kα,f)
opsα[1](elements["Γ⁴"],kα,f)

d = (k+kα)\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)

h1,l2= opsFEM[2](elements["Ωᵍ"])
L2 = log10(l2)
H1 = log10(h1)
println(L2,H1)