using  ApproxOperator, LinearAlgebra, Printf, TimerOutputs, XLSX
include("input.jl")

ndiv= 3
ndiv_p= 3
# elements,nodes,nodes_p= import_quad_GI1("./msh/square_quad_"*string(ndiv)*".msh","./msh/square_quad_"*string(ndiv_p)*".msh")
# elements,nodes,nodes_p= import_quad_GI1("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
# elements,nodes,nodes_p= import_fem_tri3_GI1("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
elements,nodes,nodes_p = import_quad_GI1_plate_with_hole("./msh/plate_with_hole_quad_"*string(ndiv)*".msh","./msh/plate_with_hole_quad_"*string(ndiv_p)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵛ"])
set∇𝝭!(elements["Ωᵛ"])
set𝝭!(elements["Γᵗ₁"])
set𝝭!(elements["Γᵗ₂"])
set𝝭!(elements["Γᵗ₃"])
set𝝭!(elements["Γᵍ₁"])
set𝝭!(elements["Γᵍ₂"])

    Ē = 3e6
    ν̄ = 0.4999999
    # ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
   
  

  T = 1000
  a = 1.0

  r(x,y) = (x^2+y^2)^0.5
  θ(x,y) = atan(y/x)
  σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
  σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
  σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))
  ApproxOperator.prescribe!(elements["Γᵗ₁"],:t₁=>(x,y,z)->σ₁₁(x,y))
  ApproxOperator.prescribe!(elements["Γᵗ₁"],:t₂=>(x,y,z)->σ₁₂(x,y))
  ApproxOperator.prescribe!(elements["Γᵗ₂"],:t₁=>(x,y,z)->σ₁₂(x,y))
  ApproxOperator.prescribe!(elements["Γᵗ₂"],:t₂=>(x,y,z)->σ₂₂(x,y))
  ApproxOperator.prescribe!(elements["Γᵗ₃"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
  ApproxOperator.prescribe!(elements["Γᵗ₃"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
  ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₂₂=>(x,y,z)->1.0)
  ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₁₁=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₁"],:n₁₂=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₁"],:g₁=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₁"],:g₂=>(x,y,z)->0.0)
  
  ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₁₁=>(x,y,z)->1.0)
  ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₂₂=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₂"],:n₁₂=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₂"],:g₁=>(x,y,z)->0.0)
  ApproxOperator.prescribe!(elements["Γᵍ₂"],:g₂=>(x,y,z)->0.0)

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

# ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
# ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
# ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
# ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
# ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
# ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
# ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->E/(1+ν)/(1-2ν)*((1-ν)*∂u∂x(x,y) + ν*∂v∂y(x,y))*n₁+E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₂)
# ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->E/(1+ν)/2*(∂u∂y(x,y) + ∂v∂x(x,y))*n₁+E/(1+ν)/(1-2ν)*(ν*∂u∂x(x,y) + (1-ν)*∂v∂y(x,y))*n₂)
# ApproxOperator.prescribe!(elements["Ω"],:b₁=>(x,y,z)->-E/(1+ν)/(1-2ν)*((1-ν)*∂²u∂x²(x,y) + ν*∂²v∂x∂y(x,y)) - E/(1+ν)/2*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y)))
# ApproxOperator.prescribe!(elements["Ω"],:b₂=>(x,y,z)->-E/(1+ν)/2*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y)) - E/(1+ν)/(1-2ν)*(ν*∂²u∂x∂y(x,y) + (1-ν)*∂²v∂y²(x,y)))


ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e9*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:Hₑ_Incompressible}(:E=>Ē,:ν=>ν̄ ),

]
opsᵛ = [
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]
opsᵈ = [
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ )
]

kᵛ = zeros(2*nᵤ,2*nᵤ)
kᵈ = zeros(2*nᵤ,2*nᵤ)
kᵍ = zeros(2*nₚ,2*nₚ)
f = zeros(2*nᵤ)

opsᵛ[1](elements["Ωᵛ"],kᵛ)
opsᵈ[1](elements["Ω"],kᵈ)  
# ops[2](elements["Ω"],f)
ops[3](elements["Γᵗ₁"],f)
ops[3](elements["Γᵗ₂"],f)
ops[3](elements["Γᵗ₃"],f)
ops[4](elements["Γᵍ₁"],kᵍ,f)
ops[4](elements["Γᵍ₂"],kᵍ,f)

d = (kᵛ+kᵈ+kᵍ)\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
prescribe!(elements["Ω"],:u=>(x,y,z)->T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) ))
    prescribe!(elements["Ω"],:v=>(x,y,z)->T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) ))
    prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
    prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))

h1,l2 = ops[6](elements["Ω"])
L2 = log10(l2)
H1 = log10(h1)
println(L2,H1)
