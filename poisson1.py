# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from netgen.meshing import *
from netgen.meshing import Element0D, Element1D, Element2D, MeshPoint, FaceDescriptor
from netgen.csg import Pnt

# generate a 1D mesh
m = Mesh(dim = 1)
nel = 10
pnums = []
for i in range(nel+1):
    pnums.append(m.Add(MeshPoint(Pnt(i/nel, 0, 0))))

# add segments
for i in range(nel):
    m.Add(Element1D([pnums[i],pnums[i+1]], index=1))

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


# Specify Dirichlet boundary conditions
fes = H1(ngsmesh, order=1, dirichlet='left|right')

#fes = H1(ngsmesh, order=1, dirichlet=[1,2])

# print ("freedofs:\n", fes.FreeDofs())

u = fes.TrialFunction()  # symbolic object
v = fes.TestFunction()   # symbolic object
gfu = GridFunction(fes)  # solution 

#a = BilinearForm(fes, symmetric=True)
a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()
#print ("mat = \n", a.mat)

# Forcing function
f = CoefficientFunction(2)

lf = LinearForm(fes)    
lf += SymbolicLFI(f*v)
lf.Assemble()

#print ("rhs = \n", lf.vec)

u = GridFunction(fes)
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * lf.vec

#print ("sol =\n", u.vec)

pnts = [i/100 for i in range(101)]
pnts_vals = [ (x,u(x)) for x in pnts if ngsmesh.Contains(x)]

#matplotlib inline
import matplotlib.pyplot as plt
pnts,vals = zip(*pnts_vals)
plt.plot(pnts,vals, "-*")
plt.show()


#L2 error 
exact = x*(1-x)
print ("L2-error:", sqrt(Integrate((u-exact)*(u-exact), ngsmesh)))