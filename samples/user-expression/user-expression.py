from dolfin import *

class MyExpression(UserExpression):

    def eval(self, values, x):
        values[0] = x[0] + x[1]

    def value_shape(self):
        return ()

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'Lagrange', 1)
f = MyExpression()
res = project(f, V)

