
using BubbleMsh
# path = "./msh/plate_with_hole.msh"
# path = "./msh/cantilever.msh"
path = "./msh/cook_membrane.msh"

#cantilever
# bubblemsh(path,[24.0,0.0,0.0],[23.0,5.0,0.0],43,2.0,0.1)

#cook_membrane
bubblemsh(path,[24.0,37.0,0.0],[23.0,7.0,0.0],20,3.0,0.05)
#plate_with_hole