
using BubbleMsh
path = "./msh/plate_with_hole.msh"
# path = "./msh/cantilever.msh"
# path = "./msh/cantilever_new_bubble.msh"
# path = "./msh/cook_membrane.msh"
# path = "./msh/patchtest.msh"
#  path = "./msh/square.msh"

 #square
#  bubblemsh(path,[0.5,0.5,0.0],[0.4,0.4,0.0],225,0.05,0.01)


#cantilever
# 4165:3845,0.225,0.1
# 3100:0.26
# 500:0.65
# 93:0.9
# bubblemsh(path,[24.0,0.0,0.0],[23.0,5.0,0.0],193,0.9,0.07)
# bubblemsh(path,[28.0,0.0,0.0],[19.5,5.5,0.0],176,0.93,0.07)
# bubblemsh(path,[4.0,0.0,0.0],[3.5,5.5,0.0],10,1.5,0.07)

#cook_membrane
# bubblemsh(path,[25.0,30.0,0.0],[20.0,9.0,0.0],55,2.9,0.1)
# plate_with_hole
# 1017,0.11

bubblemsh(path,[2.6,2.6,0.0],[1.9,1.9,0.0], 2376, 0.08, 0.05)
#square
# 127:10,80,0.6
# 420：,344,0.0285,0.05
# 54:26,0.13,0.01
# 3881,0.01,0.01
# bubblemsh(path,[0.5,0.5,0.0],[0.48,0.48,0.0],1510,0.013,0.1)