from pylab import *

def solvePDE():
    Nr, Nz = 100, 100
    Rmax = 3
    R = np.linspace(0, Rmax, Nr)

    grid=np.zeros((Nr, Nz),dtype='complex')
    grid[0, :] = np.exp(-np.square(R))

    x = 0
    while should_run(x):
        d1 = np.gradient(grid)
        dEdr = d1[1]
        dEdz = d1[0]
        dEdrr = np.gradient(dEdr)[1]
        E3 = grid**3
        q = 1+1j;
        x += 1;



    a = np.gradient(grid)
    print(a[0])
    print(a[1])


    plot(R, grid[0, :])
    plot(R, a[1][0, :])
    show()

def should_run(x):
    return x < 10

solvePDE()

