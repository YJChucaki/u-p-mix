
using ApproxOperator, LinearAlgebra, DynamicPolynomials, Random

include("import_spurious_pressure_mode.jl")

elements, nodes = import_test("msh/patchtest_11.msh")

set∇𝝭!(elements["Ω"])

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

nᵤ = 36
order = 20
@polyvar x̄ ȳ
𝒑̄ = monomials([x̄,ȳ],0:order)

∂p∂x = zeros(nᵤ)
∂p∂y = zeros(nᵤ)

# k₁: Vₙ\ker𝒫ℐₕ
# k₂: Vₙ\ker𝒫
# k₃: Vₙ\kerℐₕ𝒫
# k₄: Vₙ\ker𝒫ₕ = Vₙ\kerℐₕ-projection𝒫

k₁ = zeros(2*nᵤ,2*nᵤ)
k₂ = zeros(2*nᵤ,2*nᵤ)
k₃ = zeros(2*nᵤ,2*nᵤ)

n₁ = 100
n₂ = 100

iter = 0
while n₁ > 28 && iter <1000
global iter += 1
𝒑̃ = shuffle(𝒑̄)
∂𝒑̃∂x = differentiate.(𝒑̃[1:nᵤ],x̄)
∂𝒑̃∂y = differentiate.(𝒑̃[1:nᵤ],ȳ)
𝒑(x,y) = subs(𝒑̃[1:nᵤ], x̄=>x, ȳ=>y)
∂𝒑∂x(x,y) = subs(∂𝒑̃∂x[1:nᵤ], x̄=>x, ȳ=>y)
∂𝒑∂y(x,y) = subs(∂𝒑̃∂y[1:nᵤ], x̄=>x, ȳ=>y)

fill!(k₁,0.)
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

        # ∂p∂x .= ∂𝒑∂x(x,y)
        # ∂p∂y .= ∂𝒑∂y(x,y)
        # for (i,(∂pᵢ∂x,∂pᵢ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
        #     for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
        #         k₂[2*i-1,2*j-1] += α*∂pᵢ∂x*∂pⱼ∂x*𝑤
        #         k₂[2*i-1,2*j]   += α*∂pᵢ∂x*∂pⱼ∂y*𝑤
        #         k₂[2*i,2*j-1]   += α*∂pᵢ∂y*∂pⱼ∂x*𝑤
        #         k₂[2*i,2*j]     += α*∂pᵢ∂y*∂pⱼ∂y*𝑤
        #     end
        # end

        # fill!(∂p∂x,0.0)
        # fill!(∂p∂y,0.0)
        # N = ξ[:𝝭]
        # for (i,xᵢ) in enumerate(𝓒)
        #     ∂p∂x .+= N[i].*∂𝒑∂x(xᵢ.x,xᵢ.y)
        #     ∂p∂y .+= N[i].*∂𝒑∂y(xᵢ.x,xᵢ.y)
        # end
        # for (i,(∂pᵢ∂x,∂pᵢ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
        #     for (j,(∂pⱼ∂x,∂pⱼ∂y)) in enumerate(zip(∂p∂x,∂p∂y))
        #         k₃[2*i-1,2*j-1] += α*∂pᵢ∂x*∂pⱼ∂x*𝑤
        #         k₃[2*i-1,2*j]   += α*∂pᵢ∂x*∂pⱼ∂y*𝑤
        #         k₃[2*i,2*j-1]   += α*∂pᵢ∂y*∂pⱼ∂x*𝑤
        #         k₃[2*i,2*j]     += α*∂pᵢ∂y*∂pⱼ∂y*𝑤
        #     end
        # end
    end
end

global n₁ = rank(k₁)
# global n₂ = rank(k₂)
println(n₁)
end
alert()
# n₃ = rank(k₃)
# λ₁ = eigvals(k₁)
# λ₂ = eigvals(k₂)
# λ₃ = eigvals(k₃)