from fenics import*
import dlroms.fespaces as fe
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

mesh = fe.unitsquaremesh(30, 30)
Vh = fe.space(mesh, 'CG', 1) # FE space for u (chemical concentration)
bspace = fe.space(mesh, "CG", 1, scalar = False, bubble = True) # FE space for b (transport field)

def Stokes_solver(w1, w2):

    # Defining the mixed dimensional space for velocity and pressure
    
    pP1  = FiniteElement("CG", mesh.ufl_cell(), 1)
    vP1B = VectorElement(NodalEnrichedElement(FiniteElement("CG",       mesh.ufl_cell(), 1),
                                              FiniteElement("Bubble",   mesh.ufl_cell(), mesh.topology().dim() + 1)))

    pspace, vspace = pP1, vP1B
    W = FunctionSpace(mesh, vspace * pspace)
    (b, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)

    space = fe.space(mesh, "CG", 1, scalar = False, bubble = True)
    x, y = fe.coordinates(space).T
    f1 = -x  + 0.3
    f1[::2] = y[::2] - 0.3
    f1 = f1/((x-0.3)**2 + (y-0.3)**2+1e-8)

    f2 = -x  + 0.7
    f2[::2] = y[::2] - 0.7
    f2 = f2/((x-0.7)**2 + (y-0.7)**2+1e-8)

    f = fe.asvector(w1*f1 - w2*f2, space)

    nu = 1.0

    a = nu*inner(grad(b), grad(v))*dx - div(v)*p*dx - q*div(b)*dx
    L = inner(f, v)*dx
    dbc = DirichletBC(W.sub(0), Constant((0.0, 0.0)), lambda x, on: on)

    A = assemble(a)
    F = assemble(L)

    dbc.apply(A)
    dbc.apply(F)

    A = csr_matrix(A.array())
    F = F[:]

    bp = spsolve(A, F)
    b = bp[W.sub(0).dofmap().dofs()]    
    return b

def Convdiff_solver(w1, w2, steps = 50, dt = 5e-4):
    
    # Assembling relevant operators
    v1, v2 = TrialFunction(Vh), TestFunction(Vh)
    M = csr_matrix(assemble(v1*v2*dx).array()) # mass matrix
    S = csr_matrix(assemble(inner(grad(v1), grad(v2))*dx).array()) # stiffness (diffusion) matrix

    b = Stokes_solver(w1, w2)
    bf = fe.asvector(b, bspace)
    B = csr_matrix(assemble(inner(bf, grad(v1))*v2*dx).array()) # transport matrix
    
    # Time-stepping scheme
    def FOMstep(u0, dt, b):
        A = M + dt*S + 1000*dt*B
        F = M @ u0
        return spsolve(A, F)    


    # Initial condition
    x, y = fe.coordinates(Vh).T
    u0 = np.exp(-16*(x-0.5)**2 -16*(y-0.5)**2)
    u = [u0] # list of states in time

    # Time loop
    for n in range(steps):
        uold = u[-1]
        unew = FOMstep(uold, dt, b)
        u.append(unew)

    return np.stack(u)
