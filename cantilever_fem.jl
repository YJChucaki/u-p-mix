using  ApproxOperator, LinearAlgebra, Printf, TimerOutputs, XLSX
include("import_prescrible_ops.jl")
include("import_cantilever.jl")

# for i in 40:50
    ndiv= 4
  
    # elements,nodes= import_cantilever_fem("./msh/cantilever_tri6_"*string(ndiv)*".msh")
     elements,nodes= import_cantilever_fem("./msh/cantilever_"*string(ndiv)*".msh")
    # elements,nodes= import_cantilever_fem("./msh/cantilever_quad_"*string(ndiv)*".msh")
    # elements,nodes= import_cantilever_fem("./msh/cantilever_quad8_"*string(ndiv)*".msh")
    P = 1000
    Ē = 3e6
    # ν̄ =  0.4999999
    ν̄ = 0.3
    E = Ē/(1.0-ν̄^2)
    ν = ν̄/(1.0-ν̄)
    L = 48
    D = 12
    I = D^3/12
    EI = E*I
    
    nₘ=21
    nᵤ = length(nodes)
    eval(prescribeForGauss)
    eval(prescribeForPenalty)
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])


    ops = [
        Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
        Operator{:∫∫vᵢbᵢdxdy}(),
        Operator{:∫vᵢtᵢds}(),
        Operator{:∫vᵢgᵢds}(:α=>1e9*E),
        Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
        Operator{:Hₑ_Incompressible}(:E=>Ē,:ν=>ν̄),
        
    ]

    k = zeros(2*nᵤ,2*nᵤ)
    kᵍ = zeros(2*nᵤ,2*nᵤ)
    f = zeros(2*nᵤ)
    f = zeros(2*nᵤ)
    ops[1](elements["Ω"],k)
    ops[3](elements["Γᵗ"],f)
    ops[4](elements["Γᵍ"],kᵍ,f)

    d = (k+kᵍ)\f
    d₁ = d[1:2:2*nᵤ]
    d₂ = d[2:2:2*nᵤ]

    push!(nodes,:d₁=>d₁,:d₂=>d₂)

    h1,l2 = ops[6](elements["Ωᵍ"])
    L2 = log10(l2)
    H1 = log10(h1)
    h = log10(12.0/ndiv)
    println(L2,H1)
    # h = log10(10.0/ndiv)

#     index = 40:50
#     XLSX.openxlsx("./xlsx/mix.xlsx", mode="rw") do xf
#         Sheet = xf[2]
#         ind = findfirst(n->n==ndiv,index)+1
#         Sheet["F"*string(ind)] = h
#         Sheet["G"*string(ind)] = L2
#         Sheet["H"*string(ind)] = H1
#     end
# end
