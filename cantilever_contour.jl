using ApproxOperator, GLMakie, CairoMakie, TimerOutputs
import Gmsh: gmsh
include("input.jl")
i=72

ndiv=4
# ndiv_p=8
# elements,nodes,nodes_p = import_quad("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
elements,nodes,nodes_p,sp,xᵖ,yᵖ,zᵖ = import_fem_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements,nodes,nodes_p = import_mf_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements,nodes,nodes_p = import_fem_tri3("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv_p)*".msh")
# elements,nodes,nodes_p = import_quad("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_quad_"*string(ndiv_p)*".msh")
nᵤ = length(nodes)
nₚ = length(nodes_p)

# s =1.5*12/ndiv_p*ones(nₚ)
# 
# push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    # set𝝭!(elements["Ωᵍ"])
    # set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])

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
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P/2/I*(D^2/4-y^2))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->-P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)))
    ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2))
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
    ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

    ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
    Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e9*E),
    Operator{:Locking_ratio_mix}(:E=>Ē,:ν=>ν̄),
    Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄),
    Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
    ]
    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    kᵤₚ = zeros(2*nᵤ,nₚ)
    kₚₚ = zeros(nₚ,nₚ)
    f = zeros(2*nᵤ)

    ops[3](elements["Ω"],kᵤᵤ)
    ops[4](elements["Ω"],elements["Ωᵖ"],kᵤₚ)
    ops[5](elements["Ωᵖ"],kₚₚ)
    ops[7](elements["Γᵍ"],kᵤᵤ,f)
    ops[6](elements["Γᵗ"],f)

    k = [kᵤᵤ kᵤₚ;kᵤₚ' kₚₚ]
    f = [f;zeros(nₚ)]

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    q  = d[2*nᵤ+1:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)

𝗠 = zeros(21)
ind = 10
xs = zeros(ind)
ys = zeros(ind)
zs = zeros(ind)
color = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0, L/2, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, D/2, ind))
        indices = sp(ξ¹,ξ²,0.0)
        Nᵖ = zeros(length(indices))
        data = Dict([:xᵖ=>(1,[ξ¹]),:yᵖ=>(1,[ξ²]),:zᵖ=>(1,[0.0]),:𝝭=>(4,Nᵖ),:𝗠=>(0,𝗠)])
        𝓒 = [nodes_p[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
         p= 0.0
         for (i,xᵢ) in enumerate(𝓒)
            p  += Nᵖ[i]*xᵢ.q
        end 
        color[I,J] = p
    end
end
fig = Figure()

ax = Axis3(fig[1, 1])

hidespines!(ax)
hidedecorations!(ax)
lines!([Point(0.0, -6.0, 0.0), Point(48.0,-6.0,0.0), Point(48.0, 6.0, 0.0), Point(0, 6.0, 0.0)],color=:black)
s = surface!(ax,x,y, color=color, colormap=:coolwarm)
Colorbar(fig[2, 1], s, vertical = false)
fig