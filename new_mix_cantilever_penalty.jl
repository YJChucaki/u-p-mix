using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf,Pardiso

ndiv=8
i=200
# ndiv_p=4 
include("import_prescrible_ops.jl")                       
include("import_cantilever.jl")
include("wirteVTK.jl")
elements, nodes ,nodes_p,Ω,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cantilever_mix_internal("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*"_internal.msh","./msh/cantilever_"*string(ndiv)*"_internal.msh")
# elements, nodes ,nodes_p ,Ω= import_cantilever_mix("./msh/cantilever_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cantilever_mix("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cantilever_mix("./msh/cantilever_HR_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p = import_cantilever_T6P3("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_Q4P1("./msh/cantilever_quad_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_Q4R1("./msh/cantilever_quad_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_bubble_"*string(i)*".msh")
# elements, nodes  = import_cantilever_Q8P3("./msh/cantilever_quad8_"*string(ndiv)*".msh")
# elements, nodes  = import_cantilever_fem("./msh/cantilever_quad8_"*string(ndiv)*".msh")
    nᵤ = length(nodes)
    nₚ = length(nodes_p)
    nₑ = length(elements["Ω"])
    nₑₚ = length(Ω)
    ##for Q4P1 
    # nₚ = length(elements["Ωᵖ"])
    ##for Q8P3
    # nₚ = 3*length(elements["Ωᵖ"])
    P = 1000
    Ē = 3e6
    # Ē = 1.0
    ν̄ = 0.4999999
    # ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 48
    D = 12
    I = D^3/12
    EI = E*I
    K=Ē/3/(1-2ν̄ )
    eval(prescribeForGauss)
    eval(prescribeForPenalty)
    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
    #  set𝝭!(elements["Γᵍᵖ"])
   
    eval(opsupmix)
    kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
    kₚᵤ = zeros(nₚ,2*nᵤ)
    kₚₚ = zeros(nₚ,nₚ)
    f = zeros(2*nᵤ)
    fp= zeros(nₚ)
    opsup[3](elements["Ω"],kᵤᵤ)
    opsup[4](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
    opsup[5](elements["Ωᵖ"],kₚₚ)
    opsup[6](elements["Γᵗ"],f)

    αᵥ = 1e9
    eval(opsPenalty)
    opsα[1](elements["Γᵍ"],kᵤᵤ,f)



    k = [kᵤᵤ kₚᵤ';kₚᵤ kₚₚ]
    f = [f;fp]
    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    q  = d[2*nᵤ+1:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)


# # # exact solution contour
# K=Ē/3/(1-2ν̄ )
# G=Ē/2/(1+ν̄ )
# p̄ = zeros(nₚ)
# i=0.0
# for p in nodes_p
#     i= p.𝐼
#     ξ¹ = p.x
#     ξ² = p.y
#     ∂ū₁∂x = -P/EI*(L-ξ¹)*ξ²
#     ∂ū₁∂y = -P/6/EI*((6*L-3*ξ¹)*ξ¹ + (2+ν )*(3*ξ²^2-D^2/4))
#     ∂ū₂∂x = P/6/EI*((6*L-3*ξ¹)*ξ¹ - 3*ν *ξ²^2 + (4+5*ν )*D^2/4)
#     ∂ū₂∂y = P/EI*(L-ξ¹)*ξ²*ν 
#     ε̄₁₁ = ∂ū₁∂x
#     ε̄₂₂ = ∂ū₂∂y
#     p̄[i]= K*(ε̄₁₁+ε̄₂₂)
# end
# push!(nodes_p,:q=>p̄)

    h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = opsup[9](elements["Ωᵍ"])
    L2 = log10(l2)
    H1 = log10(h1)
    # # H1_dil = log10(h1_dil)
    # # H1_dev = log10(h1_dev)
    println(L2,H1)
    # println(H1_dil,H1_dev)
   
    eval(VTK_mix_pressure)
    # eval(VTK_mix_pressure_exact_solution)
    # eval(VTK_mix_displacement)
    # eval(VTK_Q4P1_displacement_pressure)
    # eval(VTK_T6P3_pressure)
