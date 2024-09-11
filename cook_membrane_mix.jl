
using ApproxOperator, Tensors, JLD,LinearAlgebra, GLMakie, CairoMakie, Printf
ndiv= 5
# ndiv_p=9
i=2

include("import_prescrible_ops.jl")
include("import_cook_membrane.jl")
include("wirteVTK.jl")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_tri6_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
elements, nodes ,nodes_p,Ω = import_cook_membrane_mix("./msh/cook_membrane_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(i)*".msh")
# elements, nodes = import_cook_membrane_Q4P1("./msh/cook_membrane_quad_"*string(ndiv)*".msh")
# elements, nodes ,nodes_p = import_cook_membrane_T6P3("./msh/cook_membrane_tri6_"*string(ndiv)*".msh","./msh/cook_membrane_"*string(ndiv)*".msh")
# elements, nodes = import_cook_membrane_Q8P3("./msh/cook_membrane_quad8_"*string(ndiv)*".msh")

κ = 400942
μ = 80.1938
E = 9*κ*μ/(3*κ+μ)
ν = (3*κ-2*μ)/2/(3*κ+μ)

# ν =0.499999999
# E = 70.0
# κ=E/3/(1-2ν)
# μ=E/2/(1+ν)
λ = E*ν*(1+ν)*(1-2*ν)
# ν = 0.3333
Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
Cᵢⱼᵢⱼ = E/(1+ν)/2
nᵤ = length(nodes)
nₚ = length(nodes_p)
nₑ = length(elements["Ω"])

eval(prescribeForPenalty)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])


ops = [
    Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    # # Operator{:Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    # Operator{:∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean}(:E=>E,:ν=>ν),
    # Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean_modified}(:E=>E,:ν=>ν),
    # Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean_modified}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢuᵢds}(:α=>1e9*E),
    Operator{:∫∫pCdxy_incompressible}(),
    Operator{:∫∫pJdxdy}(),
    # Operator{:∫∫p∇vxdy}(),
    Operator{:∫∫qpdxdy}(:E=>E,:ν=>ν),
]

kᵤᵤ = zeros(2*nᵤ,2*nᵤ)
kα = zeros(2*nᵤ,2*nᵤ)
kₚᵤ = zeros(nₚ,2*nᵤ)
kₚₚ = zeros(nₚ,nₚ)
f = zeros(2*nᵤ)
fp= zeros(nₚ)
fα = zeros(2*nᵤ)
fint = zeros(2*nᵤ)
fext = zeros(2*nᵤ)
d = zeros(2*nᵤ+nₚ)
Δd= zeros(2*nᵤ+nₚ)
d₁ = zeros(nᵤ)
d₂ = zeros(nᵤ)
q = zeros(nₚ)
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
        fill!(fp,0.0)
        ops[2](elements["Ω"],fint)
        f .= fext-fint

        fill!(kᵤᵤ,0.0)
        fill!(kₚᵤ,0.0)
        fill!(kₚₚ,0.0)
        
        fill!(kα,0.0)
        fill!(fα,0.0)
        ops[1](elements["Ω"],kᵤᵤ)
        ops[6](elements["Ω"],elements["Ωᵖ"],kₚᵤ)
        # ops[5](elements["Ω"],elements["Ωᵖ"],fp)
        # ops[7](elements["Ωᵖ"],kₚₚ)
        ops[4](elements["Γᵍ"],kα,fα)

        # if iter == 1
        #     Δd .= k⁻¹*(f+fα)
        # else
        #     Δd .= k⁻¹*f
        # end
        # k = [kᵤᵤ+kα kₚᵤ';λ/κ*kₚᵤ 1.0/κ*kₚₚ]
        
        k = [kᵤᵤ+kα kₚᵤ';kₚᵤ kₚₚ]
        R = [f+fα;fp]
        Δd .= k\R
        d .+= Δd 
        d₁ .+= d[1:2:2*nᵤ]
        d₂ .+= d[2:2:2*nᵤ]
        q  .+= d[2*nᵤ+1:end]
        Δdnorm = LinearAlgebra.norm(Δd)
        dnorm += Δdnorm
        err_Δd = Δdnorm/dnorm
        

        # @printf "iter = %i, err_Δf = %e, err_Δd = %e \n" iter err_Δf err_Δd
        @printf "iter = %i, err_Δd = %e \n" iter err_Δd
    end
end 

u=d₂[3]
println(u)

