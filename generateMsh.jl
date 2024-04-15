
using BenchmarkExample
import Gmsh: gmsh

n = 11
filename = "patchtest_"
BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = norder = 1,order =2)
# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true)
# gmsh.initialize()
# gmsh.open("./msh/"*filename*".msh")
# entities = getPhysicalGroups()
# nodes = get𝑿ᵢ()
