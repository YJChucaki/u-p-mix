
using BubbleMsh
# path = "./msh/plate_with_hole_new.msh"
path = "./msh/cantilever.msh"
# path = "./msh/cook_membrane.msh"

#cantilever
# 4165:3845,0.225,0.1
# 3100:0.26
# 500:0.65
# 93:0.9

bubblemsh(path,[24.0,0.0,0.0],[23.0,5.0,0.0],33,2.0,0.1)


#cook_membrane
# bubblemsh(path,[24.0,37.0,0.0],[23.0,10.0,0.0],70,2.2,0.1)
#plate_with_hole
# bubblemsh(path,[2.5,2.5,0.0],[1.8,1.8,0.0],1274,0.11,0.07)