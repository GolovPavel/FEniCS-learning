
# coding: utf-8

# In[1]:


from fenics import *
from mshr import *

import time
import matplotlib.pyplot as plt

get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


# Equation parameters
T = 2.0
num_steps = 50
dt = T / num_steps


# In[3]:


# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)

V = FunctionSpace(mesh, 'P', 1)


# In[4]:


# Define boundary condition
def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, Constant(0), boundary)


# In[5]:


# Define initial value
u_0 = Expression('exp(-a*pow(x[0], 2) - a*pow(x[1], 2))', degree=2, a=5)
u_n = interpolate(u_0, V)


# In[6]:


# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)


# In[8]:


# Create VTK file for saving solution
vtkfile = File('./solution.pvd')

u = Function(V)
t = 0
for n in range(num_steps):
    t += dt
    
    solve(a == L, u, bc)
    
    # Save to file and plot solution
    vtkfile << (u, t)
    plot(u)
    
    u_n.assign(u)

