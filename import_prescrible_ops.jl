
const lobatto3 = ([-1.0,0.0,0.0,
                    0.0,0.0,0.0,
                    1.0,0.0,0.0],[1/3,4/3,1/3])

const lobatto5 = ([-1.0,0.0,0.0,
                   -(3/7)^0.5,0.0,0.0,
                    0.0,0.0,0.0,
                    (3/7)^0.5,0.0,0.0,
                    1.0,0.0,0.0],[1/10,49/90,32/45,49/90,1/10])

const lobatto7 = ([-1.0,0.0,0.0,
                   -(5/11+2/11*(5/3)^0.5)^0.5,0.0,0.0,
                   -(5/11-2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    0.0,0.0,0.0,
                    (5/11-2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    (5/11+2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    1.0,0.0,0.0],
                  [1/21,
                   (124-7*15^0.5)/350,
                   (124+7*15^0.5)/350,
                   256/525,
                   (124+7*15^0.5)/350,
                   (124-7*15^0.5)/350,
                   1/21])

const trilobatto3 = ([0.0000000000000000,0.5000000000000000,0.0,
                      0.5000000000000000,0.0000000000000000,0.0,
                      0.5000000000000000,0.5000000000000000,0.0],
                   0.5*[1/3,1/3,1/3])

const trilobatto13 = ([1.0000000000000000,0.0000000000000000,0.0,
                       0.0000000000000000,1.0000000000000000,0.0,
                       0.0000000000000000,0.0000000000000000,0.0,
                       0.0000000000000000,0.5000000000000000,0.0,
                       0.5000000000000000,0.0000000000000000,0.0,
                       0.5000000000000000,0.5000000000000000,0.0,
                       0.0000000000000000,0.8273268353539885,0.0,
                       0.0000000000000000,0.1726731646460114,0.0,
                       0.1726731646460114,0.0000000000000000,0.0,
                       0.8273268353539885,0.0000000000000000,0.0,
                       0.8273268353539885,0.1726731646460114,0.0,
                       0.1726731646460114,0.8273268353539885,0.0,
                       0.3333333333333333,0.3333333333333333,0.0],
                   0.5*[-0.0277777777777778,
                        -0.0277777777777778,
                        -0.0277777777777778,
                         0.0296296296296297,
                         0.0296296296296297,
                         0.0296296296296297,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.4500000000000000])

opsGauss = quote
    op = Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h)
end

opsupmix = quote
    opsup = [
        Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
        Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
        Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy}(:E=>Ē,:ν=>ν̄ ),
        Operator{:∫∫p∇vdxdy}(),
        Operator{:∫∫qpdxdy}(:E=>Ē,:ν=>ν̄),
        Operator{:∫vᵢtᵢds}(),
        # Operator{:∫vᵢgᵢds}(:α=>αᵥ*E),
        Operator{:Locking_ratio_mix}(:E=>Ē,:ν=>ν̄),
        Operator{:Hₑ_up_mix}(:E=>Ē,:ν=>ν̄),
        Operator{:Hₑ_Incompressible}(:E=>E,:ν=>ν),
    ]
end

opsPenalty = quote
    opsα = [
        Operator{:∫vᵢgᵢds}(:α=>αᵥ*E),
      
    ]
end

opsNitsche = quote
    opsv = [
        Operator{:∫𝒏𝑵𝒈dΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫∇𝑴𝒏𝒂₃𝒈dΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫MₙₙθₙdΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
        Operator{:ΔMₙₛ𝒂₃𝒈_Nitsche}(:E=>E,:ν=>ν,:h=>h),
    ]
end

opsHR = quote
    opsh = [
        Operator{:∫𝒏𝑵𝒈dΓ_HR}(),
        Operator{:∫∇𝑴𝒏𝒂₃𝒈dΓ_HR}(),
        Operator{:∫MₙₙθₙdΓ_HR}(),
        Operator{:ΔMₙₛ𝒂₃𝒈_HR}(),
    ]
end