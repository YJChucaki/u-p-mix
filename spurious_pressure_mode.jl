
using ApproxOperator, LinearAlgebra, DynamicPolynomials

include("import_spurious_pressure_mode.jl")

# elements, nodes = import_test("msh/patchtest_11.msh")
elements, nodes, nodes_p = import_test_2("msh/patchtest_6.msh","msh/patchtest_bubble_28.msh")

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωₚ"])

# Linear
# 𝒑(x,y)    = [1.0,  x,  y]
# ∂𝒑∂x(x,y) = [0.0,1.0,0.0]
# ∂𝒑∂y(x,y) = [0.0,0.0,1.0]

# Cubic
# 𝒑(x,y)    = [1.0,  x,  y,x^2,x*y,y^2,  x^3,x^2*y,x*y^2,  y^3]
# ∂𝒑∂x(x,y) = [0.0,1.0,0.0,2*x,  y,0.0,3*x^2,2*x*y,  y^2,  0.0]
# ∂𝒑∂y(x,y) = [0.0,0.0,1.0,0.0,  x,2*y,  0.0,  x^2,2*x*y,3*y^2]

# Quartic
# 𝒑(x,y)    = [1.0,  x,  y,x^2,x*y,y^2,  x^3,x^2*y,x*y^2,  y^3,  x^4,  x^3*y,x^2*y^2, x*y^3, y^4]
# ∂𝒑∂x(x,y) = [0.0,1.0,0.0,2*x,  y,0.0,3*x^2,2*x*y,  y^2,  0.0,4*x^3,3*x^2*y,2*x*y^2, y^3,  0.0]
# ∂𝒑∂y(x,y) = [0.0,0.0,1.0,0.0,  x,2*y,  0.0,  x^2,2*x*y,3*y^2,  0.0,    x^3,2*x^2*y, 3*x*y^2, 4*y^3]

α = 1e7
order = 7
@polyvar x̄ ȳ
𝒑̄ = monomials([x̄,ȳ],0:order)
∂𝒑̄∂x = differentiate.(𝒑̄,x̄)
∂𝒑̄∂y = differentiate.(𝒑̄,ȳ)
𝒑(x,y) = subs(𝒑̄, x̄=>x, ȳ=>y)
∂𝒑∂x(x,y) = subs(∂𝒑̄∂x, x̄=>x, ȳ=>y)
∂𝒑∂y(x,y) = subs(∂𝒑̄∂y, x̄=>x, ȳ=>y)

nᵤ = length(𝒑(0.,0.))

∂p∂x = zeros(nᵤ)
∂p∂y = zeros(nᵤ)

# k₁: Vₙ\ker𝒫ℐₕ
# k₂: Vₙ\ker𝒫
# k₃: Vₙ\kerℐₕ𝒫
# k₄: Vₙ\ker𝒫ₕ = Vₙ\kerℐₕ-projection𝒫

k₁ = zeros(2*nᵤ,2*nᵤ)
k₂ = zeros(2*nᵤ,2*nᵤ)
k₃ = zeros(2*nᵤ,2*nᵤ)

for elm in elements["Ω"]
    𝓒 = elm.𝓒
    𝓖 = elm.𝓖
    for ξ in 𝓖
        x = ξ.x
        y = ξ.y
        𝑤 = ξ.𝑤

        fill!(∂p∂x,0.0)
        fill!(∂p∂y,0.0)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            ∂p∂x .+= B₁[i].*𝒑(xᵢ.x,xᵢ.y)
            ∂p∂y .+= B₂[i].*𝒑(xᵢ.x,xᵢ.y)
        end
        for (i,(∂pᵢ∂x,∂pᵢ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
            for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
                k₁[2*i-1,2*j-1] += α*∂pᵢ∂x*∂pⱼ∂x*𝑤
                k₁[2*i-1,2*j]   += α*∂pᵢ∂x*∂pⱼ∂y*𝑤
                k₁[2*i,2*j-1]   += α*∂pᵢ∂y*∂pⱼ∂x*𝑤
                k₁[2*i,2*j]     += α*∂pᵢ∂y*∂pⱼ∂y*𝑤
            end
        end

        ∂p∂x .= ∂𝒑∂x(x,y)
        ∂p∂y .= ∂𝒑∂y(x,y)
        for (i,(∂pᵢ∂x,∂pᵢ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
            for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
                k₂[2*i-1,2*j-1] += α*∂pᵢ∂x*∂pⱼ∂x*𝑤
                k₂[2*i-1,2*j]   += α*∂pᵢ∂x*∂pⱼ∂y*𝑤
                k₂[2*i,2*j-1]   += α*∂pᵢ∂y*∂pⱼ∂x*𝑤
                k₂[2*i,2*j]     += α*∂pᵢ∂y*∂pⱼ∂y*𝑤
            end
        end

        fill!(∂p∂x,0.0)
        fill!(∂p∂y,0.0)
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            ∂p∂x .+= N[i].*∂𝒑∂x(xᵢ.x,xᵢ.y)
            ∂p∂y .+= N[i].*∂𝒑∂y(xᵢ.x,xᵢ.y)
        end
        for (i,(∂pᵢ∂x,∂pᵢ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
            for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
                k₃[2*i-1,2*j-1] += α*∂pᵢ∂x*∂pⱼ∂x*𝑤
                k₃[2*i-1,2*j]   += α*∂pᵢ∂x*∂pⱼ∂y*𝑤
                k₃[2*i,2*j-1]   += α*∂pᵢ∂y*∂pⱼ∂x*𝑤
                k₃[2*i,2*j]     += α*∂pᵢ∂y*∂pⱼ∂y*𝑤
            end
        end
    end
end

n₁ = rank(k₁)
n₂ = rank(k₂)
n₃ = rank(k₃)


nₚ = length(nodes_p)
kₚₚ = zeros(nₚ,nₚ)
kₚᵤ = zeros(nₚ,2*nᵤ)

for elm in elements["Ωₚ"]
    𝓒 = elm.𝓒
    𝓖 = elm.𝓖
    for ξ in 𝓖
        x = ξ.x
        y = ξ.y
        𝑤 = ξ.𝑤

        N = ξ[:𝝭]
        ∂p∂x .= ∂𝒑∂x(x,y)
        ∂p∂y .= ∂𝒑∂y(x,y)

        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                kₚₚ[I,J] += N[i]*N[j]*𝑤/α
            end
            for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
                kₚᵤ[I,2*j-1] += N[i]*∂pⱼ∂x*𝑤
                kₚᵤ[I,2*j]   += N[i]*∂pⱼ∂y*𝑤
            end
        end
    end
end

k₄ = kₚᵤ'*(kₚₚ\kₚᵤ)
n₄ = rank(k₄)

λ₁ = eigvals(k₁)
λ₂ = eigvals(k₂)
λ₃ = eigvals(k₃)
λ₄ = eigvals(k₄)