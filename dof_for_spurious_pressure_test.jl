

for i in 1:100
    n = ((i+1)*(i+2)/2)^0.5
    if n/round(n) == 1
        println(i)
        println(n)
    end
end