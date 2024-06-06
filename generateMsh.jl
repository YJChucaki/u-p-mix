
using BenchmarkExample
import Gmsh: gmsh

n = 3
# filename = "patchtest_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)

# filename = "patchtest_tri6_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, order = 2)
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = norder = 1,order =2)
# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true)
# filename = "patchtest_quad8_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true, order = 2)

# filename = "cantilever_HR_"
# BenchmarkExample.Cantilever.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)

filename = "square_"
BenchmarkExample.Square.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n)

# filename = "square_quad_"
# BenchmarkExample.Square.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true)
# filename = "square_tri6_"
# BenchmarkExample.Square.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, order = 2)
# filename = "Square_quad8_"
# BenchmarkExample.Square.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n, quad = true, order = 2)