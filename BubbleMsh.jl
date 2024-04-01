
using BubbleMsh
# path = "./msh/plate_with_hole_new.msh"
# path = "./msh/cantilever.msh"
# path = "./msh/cook_membrane.msh"
path = "./msh/patchtest.msh"

#cantilever
# 4165:3845,0.225,0.1
# 3100:0.26
# 500:0.65
# 93:0.9

bubblemsh(path,[24.0,0.0,0.0],[23.0,5.0,0.0],217,0.92,0.1)


#cook_membrane
# bubblemsh(path,[24.0,37.0,0.0],[23.0,10.0,0.0],70,2.2,0.1)
#plate_with_hole
# bubblemsh(path,[2.5,2.5,0.0],[1.8,1.8,0.0],1274,0.11,0.07)
#square
# 127:10,80,0.6
bubblemsh(path,[0.5,0.5,0.0],[0.47,0.47,0.0],96,0.06,0.1)