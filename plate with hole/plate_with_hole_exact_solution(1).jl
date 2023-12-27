
E = 3E6;ν = 0.3;T = 1000;a = 1.0;b = 5.0;

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))

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
for j in 0:inte
    for i in 0:inte
        Θ = π/2 - π/4/inte*i
        R = a + (b/cos(π/4/inte*i)-a)/inte*j
        xᵢ = R*cos(Θ)
        yᵢ = R*sin(Θ)
        x[(2*inte+1)*j+i+1] = xᵢ
        y[(2*inte+1)*j+i+1] = yᵢ
        𝜎₁₁[(2*inte+1)*j+i+1] = σ₁₁(xᵢ,yᵢ)
        𝜎₂₂[(2*inte+1)*j+i+1] = σ₂₂(xᵢ,yᵢ)
        𝜎₁₂[(2*inte+1)*j+i+1] = σ₁₂(xᵢ,yᵢ)

        Θ = π/4 - π/4/inte*i
        R = a + (b/cos(Θ)-a)/inte*j
        xᵢ = R*cos(Θ)
        yᵢ = abs(R*sin(Θ))
        x[(2*inte+1)*j+inte+i+1] = xᵢ
        y[(2*inte+1)*j+inte+i+1] = yᵢ
        𝒙 = (xᵢ,yᵢ,0.0)
        uᵢ,ε₁₁,ε₂₂,ε₁₂ = get𝝐(ap,𝒙,sp)
        𝜎₁₁[(2*inte+1)*j+inte+i+1] = σ₁₁(xᵢ,yᵢ)
        𝜎₂₂[(2*inte+1)*j+inte+i+1] = σ₂₂(xᵢ,yᵢ)
        𝜎₁₂[(2*inte+1)*j+inte+i+1] = σ₁₂(xᵢ,yᵢ)
    end
end

df = DataFrame(x=x,y=y,σ₁₁=𝜎₁₁,σ₂₂=𝜎₂₂,σ₁₂=𝜎₁₂)
XLSX.openxlsx("./xlsx/plate_with_hole.xlsx", mode="rw") do xf
    name = "exact_solution"
    name∉XLSX.sheetnames(xf) ? XLSX.addsheet!(xf,name) : nothing
    XLSX.writetable!(xf[name],df)
end
