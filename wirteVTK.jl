VTK_mix_pressure = quote
 #number of an elemnt nodes
 nₑₙ=3
 fo = open("./vtk/patchtest_tri3_mix_BB_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
#  fo = open("./vtk/cook_membrane_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
 @printf fo "# vtk DataFile Version 2.0\n"
#  @printf fo "cantilever_quad4_mix\n"
#  @printf fo "cantilever_quad8_mix\n"
#  @printf fo "cantilever_tri3_mix\n"
#  @printf fo "cantilever_tri6_mix\n"
@printf fo "patchtest_tri3_mix\n"
# @printf fo "patchtest_tri6_mix\n"
 @printf fo "ASCII\n"
 @printf fo "DATASET POLYDATA\n"
 @printf fo "POINTS %i float\n" nₚ

 for p in nodes_p
    @printf fo "%f %f %f\n" p.x p.y p.z
 end
 @printf fo "POLYGONS %i %i\n" nₑₚ (nₑₙ+1)*nₑₚ
 for ap in Ω
    𝓒 = ap.𝓒
    if nₑₙ==3
     @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==4
     @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==6
     @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==8
     @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    end
 end
 @printf fo "POINT_DATA %i\n" nₚ
 @printf fo "SCALARS P float 1\n"
 @printf fo "LOOKUP_TABLE default\n"
 for p in nodes_p
    @printf fo "%f\n" p.q 
 end
# for (j,p) in enumerate(elements["Ωᵖ"])
#     ξ,  = p.𝓖
#     N = ξ[:𝝭]
#     p = 0.0
#     for (i,xᵢ) in enumerate(p.𝓒)
#         p += N[i]*xᵢ.q
#     end
#     @printf fo "%f\n" p
# end
 close(fo)
end
VTK_mix_pressure_u = quote
    #number of an elemnt nodes
    nₑₙ=3
    # fo = open("./vtk/patchtest_tri3_mix_pressure_u_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/patchtest_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    fo = open("./vtk/cantilever_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    # fo = open("./vtk/cantilever_quad4_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_quad8_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cantilever_tri6_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
   #  fo = open("./vtk/cook_membrane_tri3_mix_pressure_"*string(ndiv)*"_"*string(i)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    # @printf fo "cantilever_quad4_mix\n"
   #  @printf fo "cantilever_quad8_mix\n"
    @printf fo "cantilever_tri3_mix\n"
   #  @printf fo "cantilever_tri6_mix\n"
#    @printf fo "patchtest_tri3_mix\n"
   # @printf fo "patchtest_tri6_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nᵤ
    
  for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
       𝓒 = ap.𝓒
       if nₑₙ==3
        @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==4
        @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==6
        @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==8
        @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       end
    end
    # @printf fo "POINT_DATA %i\n" nᵤ
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS PRESSURE float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
       for ap in elements["Ω"]
           𝓒 = ap.𝓒
           𝓖 = ap.𝓖
           ε₁₁ = 0.0
           ε₂₂ = 0.0
           ε₁₂ = 0.0
           for (i,ξ) in enumerate(𝓖)
                   B₁ = ξ[:∂𝝭∂x]
                   B₂ = ξ[:∂𝝭∂y]
                   for (j,xⱼ) in enumerate(𝓒)
                       ε₁₁ += B₁[j]*xⱼ.d₁
                       ε₂₂ += B₂[j]*xⱼ.d₂
                       ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
                   end
                  
           end
           p=K*(ε₁₁+ε₂₂)
           @printf fo "%f\n" p   
       end
    close(fo)
   end
VTK_mix_displacement = quote
  #number of an elemnt nodes
  nₑₙ=3
  fo = open("./vtk/cantilever_tri3_mix_displacement_"*string(ndiv)*".vtk","w")
  @printf fo "# vtk DataFile Version 2.0\n"
  @printf fo "cantilever_tri3_mix\n"
  @printf fo "ASCII\n"
  @printf fo "DATASET POLYDATA\n"
  @printf fo "POINTS %i float\n" nᵤ
    
  for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
  @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
  for ap in elements["Ω"]
        𝓒 = ap.𝓒
        if nₑₙ==3
         @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==4
         @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==6
         @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        elseif nₑₙ==8
         @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
        end
    end
  @printf fo "POINT_DATA %i\n" nᵤ
  @printf fo "VECTORS U float\n"
  for p in nodes
     @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
  end
  close(fo)
end
VTK_T6P3_pressure = quote
     #number of an elemnt nodes
     nₑₙ=3
     fo = open("./vtk/cantilever_T6P3_pressure_"*string(ndiv)*".vtk","w")
     @printf fo "# vtk DataFile Version 2.0\n"
     @printf fo "cantilever_tri3_mix\n"
     @printf fo "ASCII\n"
     @printf fo "DATASET POLYDATA\n"
     @printf fo "POINTS %i float\n" nₚ
     
     for p in nodes_p
         @printf fo "%f %f %f\n" p.x p.y p.z
     end
     @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
     for ap in elements["Ωᵖ"]
         𝓒 = ap.𝓒
         if nₑₙ==3
          @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==4
          @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==6
          @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==8
          @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         end
     end
     @printf fo "POINT_DATA %i\n" nₚ
     @printf fo "SCALARS P float 1\n"
     @printf fo "LOOKUP_TABLE default\n"
     for p in nodes_p
         @printf fo "%f\n" p.q 
     end
     close(fo)    
 end    
VTK_Q4P1_displacement_pressure = quote
#number of an elemnt nodes
nₑₙ=4
fo = open("./vtk/cantilever_Q4P1_"*string(ndiv)*".vtk","w")
@printf fo "# vtk DataFile Version 2.0\n"
@printf fo "cantilever_Q4P1_mix\n"
@printf fo "ASCII\n"
@printf fo "DATASET POLYDATA\n"

@printf fo "POINTS %i float\n" nᵤ
for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
end
# elment type
@printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
for ap in elements["Ω"]
    𝓒 = ap.𝓒
    if nₑₙ==3
     @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==4
     @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==6
     @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    elseif nₑₙ==8
     @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
    end
end

@printf fo "POINT_DATA %i\n" nᵤ
@printf fo "VECTORS U float\n"
for p in nodes
 @printf fo "%f %f %f\n" p.d₁ p.d₂ 0.0
end

@printf fo "CELL_DATA %i\n" nₑ
@printf fo "SCALARS P float 1\n"
@printf fo "LOOKUP_TABLE default\n"
for p in 1:nₑ
    @printf fo "%f\n"  q[p] 
end

close(fo)
end