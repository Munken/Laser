from pylab import *

def solvePDE():
    Nr, Nz = 100, 100
    Rmax = 0.5
    R = np.linspace(1E-3, Rmax, Nr)

    E=np.zeros((Nr, Nz),dtype='complex')
    E[:] = 1
    gauss = np.exp(-np.square(R))
    E[0, :] = gauss

    x = 0
    while should_run(x):
        grad = np.gradient(E)

        # First derivative
        dEdr = grad[1]
        dEdz = grad[0]
        # Boundary
        dEdr[:, 0] = 0
        E[0, :] = gauss

        # Second derivative
        dEdrr = np.gradient(dEdr)[1]


        E2 = R*E * np.conj(E)
        DE = R*(1j * dEdz + dEdrr) + dEdr





        print(x)
        x += 1



    # plot(R, grid[0, :])
    plot(R, E[1, :])
    plot(R, E[2, :])
    plot(R, E[3, :])
    show()

def should_run(x):
    return x < 80

solvePDE()

