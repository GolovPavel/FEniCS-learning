from dolfin import *
import multiprocessing


class PoissonSolver1:
    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    def run_solver(self):
        self.__run_processes()
        self.__plot_results()

    def __define_mesh(self):
        self.mesh = UnitSquareMesh(30, 30)
        self.V = FunctionSpace(self.mesh, "Lagrange", 1)

    def __define_bc(self):
        # Define boundary condition
        self.u_D = Expression('1 + x[0] * x[0] + 2 * x[1] * x[1]', degree=2)
        self.bc = DirichletBC(self.V, self.u_D, self.Boundary())

    def __define_functions(self):
        # Define variational problem
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)

        self.a = dot(grad(self.u), grad(self.v)) * dx

        self.f = [Constant(-i) for i in range(10)]
        self.L = [fun * self.v * dx for fun in self.f]

    def __run_processes(self):
        pool = multiprocessing.Pool()
        pool.map(self.solve_eq, range(10))

    def solve_eq(self, i):
        self.__define_mesh()
        self.__define_bc()
        self.__define_functions()

        u = Function(self.V)
        solve(self.a == self.L[i], u, self.bc)
        File(self.mesh.mpi_comm(), 'foo%g.xml' % i) << u

    def __plot_results(self):
        self.__define_mesh()

        for i in range(10):
            u = Function(self.V, 'foo%g.xml' % i)
            plot(u)
            import matplotlib.pyplot as plt
            plt.show()


class PoissonSolver2:
    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    def run_solver(self):
        self.__define_mesh()
        self.__define_bc()
        self.__define_functions()
        self.__run_processes()

    def __define_mesh(self):
        self.mesh = UnitSquareMesh(30, 30)
        self.V = FunctionSpace(self.mesh, "Lagrange", 1)

    def __define_bc(self):
        # Define boundary condition
        self.u_D = Expression('1 + x[0] * x[0] + 2 * x[1] * x[1]', degree=2)
        self.bc = DirichletBC(self.V, self.u_D, self.Boundary())

    def __define_functions(self):
        # Define variational problem
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)

        self.a = dot(grad(self.u), grad(self.v)) * dx

        self.f = [Constant(-i) for i in range(10)]
        self.L = [fun * self.v * dx for fun in self.f]

    def __run_processes(self):
        pool = multiprocessing.Pool()
        pool.map(self.solve_eq, range(10))

    def solve_eq(self, i):
        u = Function(self.V)
        solve(self.a == self.L[i], u, self.bc)
        File(self.mesh.mpi_comm(), 'foo%g.xml' % i) << u


class PoissonSolver3:
    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    def run_solver(self):
        self.__run_processes()
        self.__plot_results()

    def __run_processes(self):
        pool = multiprocessing.Pool()
        pool.map(self.solve_eq, range(10))

    def solve_eq(self, i):
        self.mesh = UnitSquareMesh(30, 30)
        self.V = FunctionSpace(self.mesh, "Lagrange", 1)

        self.u_D = Expression('1 + x[0] * x[0] + 2 * x[1] * x[1]', degree=2)
        self.bc = DirichletBC(self.V, self.u_D, self.Boundary())

        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)

        self.a = dot(grad(self.u), grad(self.v)) * dx

        self.f = Constant(-i)
        self.L = self.f * self.v * dx

        u = Function(self.V)
        solve(self.a == self.L, u, self.bc)
        File(self.mesh.mpi_comm(), 'foo%g.xml' % i) << u

    def __plot_results(self):
        self.mesh = UnitSquareMesh(30, 30)
        self.V = FunctionSpace(self.mesh, "Lagrange", 1)

        for i in range(10):
            u = Function(self.V, 'foo%g.xml' % i)
            plot(u)
            import matplotlib.pyplot as plt
            plt.show()


def poisson_solver4(alpha):
    mesh = UnitSquareMesh(256, 256)

    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)
    bc = DirichletBC(V, Constant(0), 'on_boundary')

    a = inner(Constant(alpha) * grad(u), grad(v)) * dx
    L = inner(Constant(1), v) * dx

    uh = Function(V)
    solve(a == L, uh, bc)
    info('alpha %g -> |uh|=%g' % (alpha, uh.vector().norm('l2')))
    File(mesh.mpi_comm(), 'foo%g.xml' % alpha) << uh


if __name__ == "__main__":
    # solver = PoissonSolver1()
    # solver.run_solver()        # Ok

    # solver = PoissonSolver2()
    # solver.run_solver()        # Error! Can't pickle dolfin objects
    #
    # solver = PoissonSolver3()
    # solver.run_solver()        # Ok
    #
    # alphas = [3, 8, 9, 10]
    # pool = multiprocessing.Pool()
    # res = pool.map(poisson_solver4, alphas)  # Ok
    # pool.close()

