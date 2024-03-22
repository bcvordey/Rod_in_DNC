####### Adding L ######
# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0
# Neumann BC on x = -L and x = 0

####  Known terms #######
k = 1     # kappa
h = 1    # distance between the rod and the obstacle
l = 5    #Length of the rod
c = 1      
rhs = -1  # right hand side

#number of elements
nel = 5

from netgen.meshing import *
from netgen.meshing import Element0D, Element1D, Element2D, MeshPoint, FaceDescriptor
from netgen.csg import Pnt

# generate a 1D mesh
m = Mesh(dim = 1)
pnums = []
pnums2 = []
for i in range(nel+1):
    pnums2.append(l*((i/nel)-1))
    pnums.append(m.Add(MeshPoint(Pnt(l*((i/nel)-1), 0, 0))))

# add segments
for i in range(nel):
    m.Add(Element1D([pnums[i],pnums[(i+1)]], index=1))

m.SetMaterial(1,'material')

# add points
m.Add (Element0D(pnums[0], index=1))
m.Add (Element0D(pnums[nel], index=2))

# set boundary condition names
m.SetBCName(0,'left')
m.SetBCName(1,'right')


import ngsolve
from ngsolve import *

ngsmesh = ngsolve.Mesh(m)
# print(ngsmesh.GetBoundaries())


fes = H1(ngsmesh, order=1)
# g = 1/2*x

#Putting g on the boundaries marked non-homogeneous Dirichlet
gfu = GridFunction(fes)
# gfu.Set(g, BND)
#Draw(gfu)


#Bilinear form
u, v = fes.TnT()
a = BilinearForm(c**2*grad(u)*grad(v)*dx + c**2*k*u*v*ds(definedon='right')).Assemble()

#Linear Form
f = LinearForm(rhs*v*dx + c**2*k*h*v*ds(definedon='right')).Assemble()



#solving
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec 
# Draw(gfu)
# Draw (-grad(gfu), ngsmesh, "Flux")


# # print ("sol =\n", u.vec)

# pnts = [i/101 for i in range(101)]
# pnts_vals = [(x,gfu(x)) for x in pnts if ngsmesh.Contains(x)]





#Error 
exact = -((rhs*x*x)/(2*c**2)) - ((rhs*l*x)/(c**2)) + h + ((rhs*l)/(c**2*k))

print ("L2-error:", sqrt(Integrate((gfu-exact)*(gfu-exact), ngsmesh)))

############# EXACT SOLN #######



print(gfu.vec)
print(pnums2)

# #matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
u_vals = gfu.vec
plt.plot(pnums2,u_vals,color ='red')
plt.show()
