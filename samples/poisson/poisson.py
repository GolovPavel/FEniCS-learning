
# coding: utf-8

# In[2]:


from dolfin import *


# In[3]:


# Create mesh and define function space
mesh = UnitSquareMesh(30, 30)
V = FunctionSpace(mesh, "Lagrange", 1)


# In[4]:


# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x, on_boundary):
    return on_boundary


# In[5]:


# Define boundary condition
u_D = Expression('1 + x[0] * x[0] + 2 * x[1] * x[1]', degree=2)
bc = DirichletBC(V, u_D, boundary)


# In[6]:


# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx


# In[7]:


# Compute solution
u = Function(V)
solve(a == L, u, bc)


# In[8]:


# Save solution in VTK format
file = File("poisson.pvd")
file << u


# In[16]:


# Plot solution
plot(u)

# In[30]:


# L^2 norm error calculation
error = errornorm(u_D, u, 'L2')
print(error)


# In[22]:


# Calculate vertex values of real and machine solutions 
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)


# In[29]:


# Calculate max error
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
print(error_max)

