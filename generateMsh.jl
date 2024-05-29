
using BenchmarkExample
import Gmsh: gmsh

n = 16
filename = "patchtest_"
BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)

# filename = "patchtest_tri6_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, order = 2)
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = norder = 1,order =2)
# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true)
# filename = "patchtest_quad8_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true, order = 2)

