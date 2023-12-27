
using ApproxOperator,XLSX,TimerOutputs, YAML

to = TimerOutput()

ndiv = 6
𝒑 = "cubic"
config = YAML.load_file("./yml/plate_with_hole_hrrk_"*𝒑*".yml")
elements, nodes = importmsh("./msh/plate_with_hole_"*string(ndiv)*".msh",config)
# elements, nodes = importmsh("./msh/plate_with_hole_irregular.msh")
nₚ = getnₚ(elements["Ω"])

set_memory_𝗠!(elements["Ω̃"],:∇̃)
set_memory_𝗠!(elements["Γᵍ₁"],:𝝭,:∇̃)
set_memory_𝗠!(elements["Γᵍ₂"],:𝝭,:∇̃)
# s = 2.5*1.5/ndiv*ones(nₚ)
# s = 3.5*1.53/ndiv*ones(nₚ)
s = 3.5*1.51/ndiv*ones(nₚ)
# s = 3.5*2/3*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
@timeit to "Total Time" begin
elements["Ω∩Γᵍ₁"] = elements["Ω"]∩elements["Γᵍ₁"]
elements["Ω∩Γᵍ₂"] = elements["Ω"]∩elements["Γᵍ₂"]
@timeit to "prescribling" begin
E = 3E6;ν = 0.3;T = 1000;a = 1.0;b = 5.0;

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))
prescribe!(elements["Γᵗ₁"],:t₁=>(x,y,z)->σ₁₁(x,y))
prescribe!(elements["Γᵗ₁"],:t₂=>(x,y,z)->σ₁₂(x,y))
prescribe!(elements["Γᵗ₂"],:t₁=>(x,y,z)->σ₁₂(x,y))
prescribe!(elements["Γᵗ₂"],:t₂=>(x,y,z)->σ₂₂(x,y))
prescribe!(elements["Γᵗ₃"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ₃"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ₁"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ₂"],:n₁₁=>(x,y,z)->1.0)
end

@timeit to "shape functions Ω" set𝝭!(elements["Ω"])
@timeit to "shape functions Ω" set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
@timeit to "shape functions Γᵍ" set∇̃𝝭!(elements["Γᵍ₁"],elements["Ω∩Γᵍ₁"])
@timeit to "shape functions Γᵍ" set∇̃𝝭!(elements["Γᵍ₂"],elements["Ω∩Γᵍ₂"])
@timeit to "shape functions Γᵍ" set∇̄𝝭!(elements["Γᵍ₁"])
@timeit to "shape functions Γᵍ" set∇̄𝝭!(elements["Γᵍ₂"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₁"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₂"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₃"])

coefficient = (:E=>E,:ν=>ν)
ops = [Operator(:∫∫εᵢⱼσᵢⱼdxdy,coefficient...),
       Operator(:∫vᵢtᵢds,coefficient...),
       Operator(:∫σᵢⱼnⱼgᵢds,coefficient...),
       Operator(:∫σ̄ᵢⱼnⱼgᵢds,coefficient...),
       Operator(:Hₑ_PlaneStress,coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

@timeit to "assembly in Ω" ops[1](elements["Ω̃"],k)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₁"],f)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₂"],f)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₃"],f)
@timeit to "assembly in Γᵍ" ops[3](elements["Γᵍ₁"],k,f)
@timeit to "assembly in Γᵍ" ops[3](elements["Γᵍ₂"],k,f)
@timeit to "assembly in Γᵍ" ops[4](elements["Γᵍ₁"],k,f)
@timeit to "assembly in Γᵍ" ops[4](elements["Γᵍ₂"],k,f)

@timeit to "solve" d = k\f
end
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set∇𝝭!(elements["Ω̄"])
prescribe!(elements["Ω̄"],:u=>(x,y,z)->T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) ))
prescribe!(elements["Ω̄"],:v=>(x,y,z)->T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) ))
prescribe!(elements["Ω̄"],:∂u∂x=>(x,y,z)->T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
prescribe!(elements["Ω̄"],:∂u∂y=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
prescribe!(elements["Ω̄"],:∂v∂x=>(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
prescribe!(elements["Ω̄"],:∂v∂y=>(x,y,z)->T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
h1,l2 = ops[5](elements["Ω"])
l2 = log10(l2)
h1 = log10(h1)
h = log10(1.0/ndiv)

show(to)

# if ndiv == 5
# Cᵢᵢᵢᵢ = E/(1-ν^2)
# Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
# Cᵢⱼᵢⱼ = E/2/(1+ν)
# inte = 100
# n̄ₚ = (inte+1)*(2*inte+1)
# x = zeros(n̄ₚ)
# y = zeros(n̄ₚ)
# 𝜎₁₁ = zeros(n̄ₚ)
# 𝜎₂₂ = zeros(n̄ₚ)
# 𝜎₁₂ = zeros(n̄ₚ)
# 𝗠 = elements["Ω"][1].𝗠
# 𝝭 = elements["Ω"][1].𝝭
# ap = ReproducingKernel{type...,:Node}([Node(i,nodes) for i in 1:nₚ],Node[],𝗠,𝝭)
# for j in 0:inte
#     for i in 0:inte
#         Θ = π/2 - π/4/inte*i
#         R = a + (b/cos(π/4/inte*i)-a)/inte*j
#         xᵢ = R*cos(Θ)
#         yᵢ = R*sin(Θ)
#         x[(2*inte+1)*j+i+1] = xᵢ
#         y[(2*inte+1)*j+i+1] = yᵢ
#         𝒙 = (xᵢ,yᵢ,0.0)
#         uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
#         𝜎₁₁[(2*inte+1)*j+i+1] = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
#         𝜎₂₂[(2*inte+1)*j+i+1] = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
#         𝜎₁₂[(2*inte+1)*j+i+1] = Cᵢⱼᵢⱼ*ε₁₂

#         Θ = π/4 - π/4/inte*i
#         R = a + (b/cos(Θ)-a)/inte*j
#         xᵢ = R*cos(Θ)
#         yᵢ = abs(R*sin(Θ))
#         x[(2*inte+1)*j+inte+i+1] = xᵢ
#         y[(2*inte+1)*j+inte+i+1] = yᵢ
#         𝒙 = (xᵢ,yᵢ,0.0)
#         uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
#         𝜎₁₁[(2*inte+1)*j+inte+i+1] = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
#         𝜎₂₂[(2*inte+1)*j+inte+i+1] = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
#         𝜎₁₂[(2*inte+1)*j+inte+i+1] = Cᵢⱼᵢⱼ*ε₁₂
#     end
# end

# df = DataFrame(σ₁₁=𝜎₁₁,σ₂₂=𝜎₂₂,σ₁₂=𝜎₁₂)
# XLSX.openxlsx("./xlsx/plate_with_hole.xlsx", mode="rw") do xf
#     name = "rigsi_hr"
#     name∉XLSX.sheetnames(xf) ? XLSX.addsheet!(xf,name) : nothing
#     XLSX.writetable!(xf[name],df)
# end
# end
