
using BenchmarkExample
import Gmsh: gmsh

n = 21
filename = "patchtest_"
BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)

# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true)
# gmsh.initialize()
# gmsh.open("./msh/"*filename*".msh")
# entities = getPhysicalGroups()
# nodes = get𝑿ᵢ()
