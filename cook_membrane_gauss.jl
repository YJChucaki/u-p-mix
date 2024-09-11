
using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf
include("input.jl")
# for i in 2:10
ndiv= 5



include("import_prescrible_ops.jl")
include("import_cook_membrane.jl")
include("wirteVTK.jl")
elements, nodes  = import_cook_membrane_fem("./msh/cook_membrane_"*string(ndiv)*".msh")
# elements, nodes  = import_cook_membrane_MF("./msh/cook_membrane_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_quad_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
nᵤ = length(nodes)
nₑ = length(elements["Ω"])

κ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
ν = (3*κ-2*μ)/2/(3*κ+μ)
# ν =0.499999999
# E = 70.0
# ν = 0.3333
Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
Cᵢⱼᵢⱼ = E/(1+ν)/2

eval(prescribeForPenalty)

set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])




ops = [
    Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢuᵢds}(:α=>1e7*E),
]

k = zeros(2*nᵤ,2*nᵤ)
kα = zeros(2*nᵤ,2*nᵤ)
f = zeros(2*nᵤ)
fα = zeros(2*nᵤ)
fint = zeros(2*nᵤ)
fext = zeros(2*nᵤ)
d = zeros(2*nᵤ)
Δd= zeros(2*nᵤ)
d₁ = zeros(nᵤ)
d₂ = zeros(nᵤ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)

nmax = 20
P = 0:6.25/nmax:6.25
tolerance=1.0e-10;maxiters=1000;
for (n,p) in enumerate(P)
    if n == 1
        continue
    end
    err_Δd = 1.0
    dnorm = 0.0
    # err_Δf = 1.0
    # fnorm = 0.0
    @printf "Load step=%i,p=%e \n" n p
    fill!(fext,0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->p)
    ops[3](elements["Γᵗ"],fext)
    # fill!(k,0.0)
    # fill!(kα,0.0)
    # fill!(fα,0.0)
    # ops[1](elements["Ω"],k)
    # ops[4](elements["Γᵍ"],kα,fα)
    # k⁻¹ .= inv(k+kα)
    iter = 0
    while err_Δd>tolerance && iter<maxiters
        iter += 1
        fill!(fint,0.0)
        ops[2](elements["Ω"],fint)
        f .= fext-fint

        fill!(k,0.0)
        fill!(kα,0.0)
        fill!(fα,0.0)
        ops[1](elements["Ω"],k)
        ops[4](elements["Γᵍ"],kα,fα)

        # if iter == 1
        #     Δd .= k⁻¹*(f+fα)
        # else
        #     Δd .= k⁻¹*f
        # end

        Δd .= (k+kα)\(f+fα)

        # fnorm = norm(f)
        # fᵗnorm = fnorm+1.0
        # fᵗ .= f
        # λ = 2.0
        # while fᵗnorm ≥ fnorm && λ > tolerance
        #     # println(λ)
        #     fill!(fint,0.0)
        #     λ *= 0.5
        #     d₁ .= d[1:2:2*nₚ]+λ*Δd[1:2:2*nₚ]
        #     d₂ .= d[2:2:2*nₚ]+λ*Δd[2:2:2*nₚ]
        #     ops[2](elements["Ω"],fint)
        #     fᵗ = fext-fint
        #     fᵗnorm = norm(fᵗ)
        #     # println(fnorm)
        #     # println(fᵗnorm)
        # end
        # d .+= λ*Δd 

        d .+= Δd 
        d₁ .= d[1:2:2*nᵤ]
        d₂ .= d[2:2:2*nᵤ]

        Δdnorm = LinearAlgebra.norm(Δd)
        # Δdnorm = LinearAlgebra.norm(λ*Δd)
        dnorm += Δdnorm
        err_Δd = Δdnorm/dnorm
        # err_Δd = Δdnorm
        # Δfnorm = LinearAlgebra.norm(f+fα)
        # fnorm += Δfnorm
        # err_Δf = Δfnorm/fnorm

        # @printf "iter = %i, err_Δf = %e, err_Δd = %e \n" iter err_Δf err_Δd
        @printf "iter = %i, err_Δd = %e \n" iter err_Δd
    end
end 

u=d₂[3]
println(u)

