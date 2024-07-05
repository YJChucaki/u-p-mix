
using ApproxOperator, Tensors,  LinearAlgebra, Printf
include("import_heat_conduction.jl")
include("wirteVTK.jl")
# for i=2:10
   
ndiv =17
# println(nₚ)
elements,nodes = import_patchtest_fem("./msh/patchtest_tri6_"*string(ndiv)*".msh")
nᵤ = length(nodes)
nₑ = length(elements["Ω"])
# nₑₚ = length(Ω)


set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])
D=1   #thermal conductivity coefficient
t=1 #thickness

n = 1
T(x,y) = (x+y)^n
∂T∂x(x,y) = n*(x+y)^abs(n-1)
∂T∂y(x,y) = n*(x+y)^abs(n-1)
∂²T∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²T∂y²(x,y) = n*(n-1)*(x+y)^abs(n-2)
s(x,y) = -D*(∂²T∂x²(x,y)+∂²T∂y²(x,y))

eval(prescribeForFem)

ops = [
       Operator{:∫Tᵢhᵢds}(:t=>t),
       Operator{:∫Tᵢgᵢds}(:α=>1e12*D,:t=>t),
       Operator{:∫∫Tᵢsᵢdxdy}(:t=>t),
       Operator{:T_error}(:D=>D,:t=>t),
       Operator{:∫∫∇TᵢD∇Tⱼdxdy}(:D=>D,:t=>t),
]


kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
k = zeros(nᵤ,nᵤ)
f = zeros(nᵤ)


ops[5](elements["Ω"],k)
ops[3](elements["Ω"],f)
ops[2](elements["Γ"],kᵅ,f)

q = (k+kᵅ)\f #temperatures

push!(nodes,:T=>q)


set𝝭!(elements["Ωᵍ"])
l2= ops[4](elements["Ωᵍ"])
# h1,l2,h1_dil,h1_dev= ops[5](elements["Ωᵍ"],elements["Ωᵍᵖ"])
L2 = log10(l2)
# H1 = log10(h1)
# H1_dil = log10(h1_dil)
# H1_dev = log10(h1_dev)
           
println(L2)

# eval(VTK_mix_pressure)


