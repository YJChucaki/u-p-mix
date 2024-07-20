using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf,Pardiso
# NP=[40,80,120,140]
# for n=1:4
    # i=NP[n]
ndiv= 6
#  ndiv_p=8
i= 20
# 40,60-3
# 80-4
# 100,120-5
# 160,200-7

include("import_prescrible_ops.jl")
include("import_cantilever.jl")
include("wirteVTK.jl")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_quad_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_quad8_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
elements, nodes ,nodes_p ,Ω,xᵖ,yᵖ,zᵖ, sp,type = import_cantilever_mix("./msh/square_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes ,nodes_p ,xᵖ,yᵖ,zᵖ, sp,type= import_cantilever_mix("./msh/square_tri6_"*string(ndiv)*".msh","./msh/patchtest_bubble_"*string(i)*".msh")
# elements, nodes = import_cantilever_Q4P1("./msh/square_quad_"*string(ndiv)*".msh")
# elements, nodes = import_cantilever_Q8P3("./msh/square_quad8_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p = import_cantilever_mix("./msh/square_tri6_"*string(ndiv)*".msh","./msh/square_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p = import_cantilever_T6P3("./msh/square_tri6_"*string(ndiv)*".msh","./msh/square_"*string(ndiv)*".msh")   

nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑ = length(elements["Ω"])
nₑₚ = length(Ω)
    ##for Q4P1 
    # nₚ = length(elements["Ωᵖ"])
    ##for Q8P3 
    # nₚ = 3*length(elements["Ωᵖ"])
    nₘ=21
    P = 1000
    Ē = 3e6
    # Ē = 1.0
    # ν̄ = 0.4999999
    ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 10
    D = 10
    I = D^3/12
    EI = E*I
    K=Ē/3/(1-2ν̄ )
    eval(prescribeForSquare)
    set𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
    # set𝝭!(elements["Γᵍᵖ"])
   

    

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
    # opsα[2](elements["Γᵍ"],elements["Γᵍᵖ"],kᵤₚ,fp)

    
    k = [kᵤᵤ kₚᵤ';kₚᵤ kₚₚ]
    f = [f;fp]

    d = k\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]
    q  = d[2*nᵤ+1:end]
    push!(nodes,:d₁=>d₁,:d₂=>d₂)
    push!(nodes_p,:q=>q)

    # kᵈ = kᵤᵤ
    # kᵛ = -kₚᵤ'*(kₚₚ\kₚᵤ)
    # vᵈ = eigvals(kᵈ)
    # vᵛ = eigvals(kᵛ)
    # γ = eigvals(kᵛ,kᵈ)
    # println(γ[2*nᵤ-nₚ+1])

    h1,l2,h1_dil,h1_dev = opsup[8](elements["Ωᵍ"],elements["Ωᵍᵖ"])
    # h1,l2 = opsup[8](elements["Ω"],elements["Ωᵖ"])
    L2 = log10(l2)
    H1 = log10(h1)
    # H1_dil = log10(h1_dil)
    # H1_dev = log10(h1_dev)
   
    println(L2,H1)
    # println(H1_dil,H1_dev)
    # println(l2,h1)
    # println(h1_dil,h1_dev)
    # h = log10(10.0/ndiv)
    eval(VTK_mix_pressure)