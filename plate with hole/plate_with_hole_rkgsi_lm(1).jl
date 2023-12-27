
using Revise,ApproxOperator,DataFrames,TimerOutputs,XLSX

to = TimerOutput()

ndiv = 6
elements, nodes = importmsh("./msh/plate_with_hole_"*string(ndiv)*".msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

# type = (SNode,:Quadratic2D,:□,:CubicSpline)
# s = 2.5*1.5/ndiv*ones(nₚ)
type = (SNode,:Cubic2D,:□,:CubicSpline)
s = 3.5*1.51/ndiv*ones(nₚ)
sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 5)
@timeit to "Total Time" begin
@timeit to "searching" begin
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],sp)
elements["Ω̃"] = ReproducingKernel{type...,:Tri3}(elements["Ω"])
elements["Γᵗ₁"] = ReproducingKernel{type...,:Seg2}(elements["Γᵗ₁"],sp)
elements["Γᵗ₂"] = ReproducingKernel{type...,:Seg2}(elements["Γᵗ₂"],sp)
elements["Γᵗ₃"] = ReproducingKernel{type...,:Seg2}(elements["Γᵗ₃"],sp)
elements["Γˡ"] = Element{:Seg2}([elements["Γᵍ₁"];elements["Γᵍ₂"]],renumbering=true)
elements["Γᵍ₁"] = ReproducingKernel{type...,:Seg2}(elements["Γᵍ₁"],sp)
elements["Γᵍ₂"] = ReproducingKernel{type...,:Seg2}(elements["Γᵍ₂"],sp)
nₗ = getnₚ(elements["Γˡ"])

end
@timeit to "prescribling" begin

# set𝓖!(elements["Ω"],:TriRK6,:∂1)
# set𝓖!(elements["Ω̃"],:TriGI3,:∂1,:∂x,:∂y,:∂z)
# set𝓖!(elements["Γᵗ₁"],:SegRK3,:∂1)
# set𝓖!(elements["Γᵗ₂"],:SegRK3,:∂1)
# set𝓖!(elements["Γᵗ₃"],:SegRK3,:∂1)
# set𝓖!(elements["Γᵍ₁"],:SegRK3,:∂1)
# set𝓖!(elements["Γᵍ₂"],:SegRK3,:∂1)
# set𝓖!(elements["Γ̃ᵍ₁"],:SegRK3,:∂1,:∂x,:∂y,:∂z)
# set𝓖!(elements["Γ̃ᵍ₂"],:SegRK3,:∂1,:∂x,:∂y,:∂z)

set𝓖!(elements["Ω"],:TriRK13,:∂1)
set𝓖!(elements["Ω̃"],:TriGI6,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["Γᵗ₁"],:SegRK5,:∂1)
set𝓖!(elements["Γᵗ₂"],:SegRK5,:∂1)
set𝓖!(elements["Γᵗ₃"],:SegRK5,:∂1)
set𝓖!(elements["Γᵍ₁"],:SegRK5,:∂1)
set𝓖!(elements["Γᵍ₂"],:SegRK5,:∂1)

E = 3E6;ν = 0.3;T = 1000;a = 1.0;

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))
prescribe!(elements["Γᵗ₁"],:t₁,(x,y,z)->σ₁₁(x,y))
prescribe!(elements["Γᵗ₁"],:t₂,(x,y,z)->σ₁₂(x,y))
prescribe!(elements["Γᵗ₂"],:t₁,(x,y,z)->σ₁₂(x,y))
prescribe!(elements["Γᵗ₂"],:t₂,(x,y,z)->σ₂₂(x,y))
prescribe!(elements["Γᵗ₃"],:t₁,(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ₃"],:t₂,(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ₁"],:n₂₂,(x,y,z)->1.0)
prescribe!(elements["Γᵍ₂"],:n₁₁,(x,y,z)->1.0)
end

push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
@timeit to "shape functions Ω" set𝝭!(elements["Ω"])
@timeit to "shape functions Ω" set∇̃𝝭!(elements["Ω̃"],elements["Ω"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₁"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₂"])
@timeit to "shape functions Γᵗ" set𝝭!(elements["Γᵗ₃"])
@timeit to "shape functions Γᵍ" set𝝭!(elements["Γᵍ₁"])
@timeit to "shape functions Γᵍ" set𝝭!(elements["Γᵍ₂"])

coefficient = (:E=>E,:ν=>ν,:α=>1E3*E)
ops = [Operator(:∫∫εᵢⱼσᵢⱼdxdy,coefficient...),
       Operator(:∫vᵢtᵢds,coefficient...),
       Operator(:∫λᵢgᵢds,coefficient...),
       Operator(:Hₑ_PlaneStress,coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
g = zeros(2*nₗ,2*nₚ)
q = zeros(2*nₗ)

@timeit to "assembly in Ω" ops[1](elements["Ω̃"],k)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₁"],f)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₂"],f)
@timeit to "assembly in Γᵗ" ops[2](elements["Γᵗ₃"],f)
@timeit to "assembly in Γᵍ" ops[3]([elements["Γᵍ₁"];elements["Γᵍ₂"]],elements["Γˡ"],g,q)
nᵍ₁ = length(elements["Γᵍ₁"])
nᵍ₂ = length(elements["Γᵍ₁"])
index = [collect(2:2:nₗ);collect(nₗ+1:2:2*nₗ-1)]

@timeit to "solve" d = [k g[index,:]';g[index,:] zeros(nₗ,nₗ)]\[f;q[index]]
end
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
set𝓖!(elements["Ω"],:TriGI16,:∂1,:∂x,:∂y,:∂z)
set∇𝝭!(elements["Ω"])
prescribe!(elements["Ω"],:u,(x,y,z)->T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) ))
prescribe!(elements["Ω"],:v,(x,y,z)->T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) ))
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
prescribe!(elements["Ω"],:∂v∂x,(x,y,z)->T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y))))
prescribe!(elements["Ω"],:∂v∂y,(x,y,z)->T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y))))
h1,l2 = ops[4](elements["Ω"])
l2 = log10(l2)
h1 = log10(h1)
h = log10(1.0/ndiv)

show(to)

#
if ndiv == 6
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)
inte = 100
n̄ₚ = (inte+1)*(2*inte+1)
x = zeros(n̄ₚ)
y = zeros(n̄ₚ)
𝜎₁₁ = zeros(n̄ₚ)
𝜎₂₂ = zeros(n̄ₚ)
𝜎₁₂ = zeros(n̄ₚ)
𝗠 = elements["Ω"][1].𝗠
𝝭 = elements["Ω"][1].𝝭
ap = ReproducingKernel{type...,:Node}([Node(i,nodes) for i in 1:nₚ],Node[],𝗠,𝝭)
for j in 0:inte
    for i in 0:inte
        Θ = π/2 - π/4/inte*i
        R = a + (b/cos(π/4/inte*i)-a)/inte*j
        xᵢ = R*cos(Θ)
        yᵢ = R*sin(Θ)
        x[(2*inte+1)*j+i+1] = xᵢ
        y[(2*inte+1)*j+i+1] = yᵢ
        𝒙 = (xᵢ,yᵢ,0.0)
        uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
        𝜎₁₁[(2*inte+1)*j+i+1] = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        𝜎₂₂[(2*inte+1)*j+i+1] = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
        𝜎₁₂[(2*inte+1)*j+i+1] = Cᵢⱼᵢⱼ*ε₁₂

        Θ = π/4 - π/4/inte*i
        R = a + (b/cos(Θ)-a)/inte*j
        xᵢ = R*cos(Θ)
        yᵢ = abs(R*sin(Θ))
        x[(2*inte+1)*j+inte+i+1] = xᵢ
        y[(2*inte+1)*j+inte+i+1] = yᵢ
        𝒙 = (xᵢ,yᵢ,0.0)
        uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
        𝜎₁₁[(2*inte+1)*j+inte+i+1] = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        𝜎₂₂[(2*inte+1)*j+inte+i+1] = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂
        𝜎₁₂[(2*inte+1)*j+inte+i+1] = Cᵢⱼᵢⱼ*ε₁₂
    end
end

df = DataFrame(σ₁₁=𝜎₁₁,σ₂₂=𝜎₂₂,σ₁₂=𝜎₁₂)
XLSX.openxlsx("./xlsx/plate_with_hole.xlsx", mode="rw") do xf
    name = "rigsi_lm"
    name∉XLSX.sheetnames(xf) ? XLSX.addsheet!(xf,name) : nothing
    XLSX.writetable!(xf[name],df)
end
end
