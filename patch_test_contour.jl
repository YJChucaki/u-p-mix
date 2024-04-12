using ApproxOperator, JLD, GLMakie, Tensors

import Gmsh: gmsh

include("import_patchtest.jl")
# for i=2:10
   
ndiv= 11
nₚ = 105
# println(nₚ)
# elements,nodes,nodes_p = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")
elements,nodes,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(nₚ)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)
 
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γ"])
Ē = 1.0
ν̄ = 0.4999999
# ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

# n = 1
# u(x,y) = (x+y)^n
# v(x,y) = (x+y)^n
# ∂u∂x(x,y) = n*(x+y)^abs(n-1)
# ∂u∂y(x,y) = n*(x+y)^abs(n-1)
# ∂v∂x(x,y) = n*(x+y)^abs(n-1)
# ∂v∂y(x,y) = n*(x+y)^abs(n-1)
# ∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²v∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
n = 2
u(x,y) = (1+2*x+3*y)^n
v(x,y) = (4+5*x+6*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))
∂σ₁₁∂x(x,y) = E/(1-ν^2)*(∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1-ν^2)*(∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂x(x,y) + ∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂y(x,y) + ∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)

eval(prescribe)

ops = [
       Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
       Operator{:∫vᵢtᵢds}(),
       Operator{:∫vᵢgᵢds}(:α=>1e13*E),
       Operator{:∫∫vᵢbᵢdxdy}(),
       Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄)
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
ops[3](elements["Γ"],kᵤᵤ,f)
ops[4](elements["Ω"],f)


# kᵈ = kᵤᵤ
# kᵛ = kᵤₚ*(kₚₚ\kᵤₚ')
k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
f = [f;zeros(nₚ)]
# d = (kᵛ+kᵈ)\f

d = k\f
d₁ = d[1:2:2*nᵤ]
d₂ = d[2:2:2*nᵤ]
p  = d[2*nᵤ+1:end]

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes_p,:q=>p)


# # exact solution contour
K=Ē/3/(1-2ν̄ )
G=Ē/2/(1+ν̄ )
𝗠 = zeros(21)
ind = 20
xs = zeros(ind)
ys = zeros(ind)
color = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0,1.0, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 1.0, ind))
        ∂ū₁∂x  = 2*n*(1+2*ξ¹+3*ξ²)^abs(n-1)
        ∂ū₁∂y = 3*n*(1+2*ξ¹+3*ξ²)^abs(n-1)
        ∂ū₂∂x = 5*n*(4+5*ξ¹+6*ξ²)^abs(n-1)
        ∂ū₂∂y = 6*n*(4+5*ξ¹+6*ξ²)^abs(n-1)
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        xs[I] = ξ¹
        ys[J] = ξ² 
        color[I,J] = K*(ε̄₁₁+ε̄₂₂)
    end
end


# 𝗠 = zeros(21)
# ind = 20
# xs = zeros(ind)
# ys = zeros(ind)
# color = zeros(ind,ind)

# for (I,ξ¹) in enumerate(LinRange(0.0, 1.0, ind))
#     for (J,ξ²) in enumerate(LinRange(0.0, 1.0, ind))
#         indices = sp(ξ¹,ξ²,0.0)
#         Nᵖ = zeros(length(indices))
#         data = Dict([:x=>(1,[ξ¹]),:y=>(1,[ξ²]),:z=>(1,[0.0]),:𝝭=>(4,Nᵖ),:𝗠=>(0,𝗠)])
#         𝓒 = [nodes_p[k] for k in indices]
#         𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
#         ap = type(𝓒,𝓖)
#         set𝝭!(ap)
#          p= 0.0       
#         for (i,xᵢ) in enumerate(𝓒)
#             p  += Nᵖ[i]*xᵢ.q
           
#         end 
#         xs[I] = ξ¹
#         ys[J] = ξ² 
#         color[I,J] = p
        
#     end
# end

fig = Figure()
ax = Axis(fig[1, 1])
hidespines!(ax)
hidedecorations!(ax)

# s=surface!(xs,ys, color, colormap=:coolwarm)
# s = contourf!(xs,ys, color, colormap=:coolwarm,levels=-1000:200:1000)
s = contourf!(xs,ys, color, colormap=:coolwarm)
Colorbar(fig[1, 2], s)

# # elements
lwb = 2.5;lwm =2.5;mso =5;msx =15;ppu = 2.5;α = 0.7;
for elm in elements["Ω"]
   
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
   
    lines!(x,y, linewidth = 0.3, color = :black)

end
# scatter!(x,y,marker = :circle, markersize = mso, color = :black)
lines!([0.0,1.0,1.0,0.0,0.0],[0.0,0.0,1.0,1.0,0.0], linewidth = lwb, color = :black)
# save("./png/patchtest_tri3_"*string(nₚ)*".png",fig)
save("./png/patchtest_tri3_exact.png",fig)
fig
# end