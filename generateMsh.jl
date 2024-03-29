
using BenchmarkExample
import Gmsh: gmsh

n = 11
filename = "patchtest"
BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*".msh", transfinite = n)

# gmsh.initialize()
# gmsh.open("./msh/"*filename*".msh")
# entities = getPhysicalGroups()
# nodes = get𝑿ᵢ()
