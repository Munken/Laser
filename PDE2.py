from pylab import *

plot3d = True;

try:
     from mpl_toolkits.mplot3d import Axes3D
except ImportError:
     plot3D = False;
     print('Not plotting 3D.')



Nr = 500; Nz = 100000
Rmax = 2.5; Zmax = 0.005
dR = Rmax / Nr; dZ = Zmax / Nz
R = np.linspace(1E-3, Rmax, Nr)
R2 = dR ** 2


def dEdr2C(E, i, n):
    return (E[n - 1, i + 1] - 2 * E[n - 1, i] + E[n - 1, i - 1]) / R2


def dEdrC(E, i, n):
    return R[i] * (E[n - 1, i + 1] - E[n - 1, i - 1] / dR)


def E3(E, i, n):
    return (np.abs(E[n - 1, i]) ** 2) * E[n - 1, i]


def dEdr2R(E, n):
    return (E[n - 1, 2] + E[n - 1, 1] - 2 * E[n - 1, 0]) / R2


def dEdr2L(E, n):
    return (E[n - 1, -3] + E[n - 1, -1] - 2 * E[n - 1, -2]) / R2




def solve_PDE2():
    E=np.zeros((Nz+1, Nr),dtype='complex')
    E[0,:] = 1E-3*np.exp(-np.square(R))


    n = 1
    for i in range(1, Nr-1):
        E[n,i] = 1j*dZ*(dEdr2C(E, i, n) + dEdrC(E, i, n) + np.abs(E[n-1,i])**2*E[n-1,i]) + E[n-1, i]

    E[n,0] = 1j*dZ*(dEdr2R(E, n) + np.abs(E[n-1,0])**2*E[n-1,0]) + E[n-1,0]
    E[n,Nr-1] = 1j*dZ*(dEdr2L(E, n) + R[-1]*(E[n-1, -2] - E[n-1, -1]/dR) + E3(E, -1, n)) + E[n-1,-1]

    for n in range(2, Nz+1):
        for i in range(1, Nr-1):
             E[n,i] = (2j*dZ*(dEdr2C(E, i, n) + dEdrC(E, i, n) + E3(E, i, n)) + 4*E[n-1, i] - E[n-2, i])/3

        E[n,0] = (2j*dZ*(dEdr2R(E, n) + np.abs(E[n-1,0])**2*E[n-1,0]) + 4*E[n-1, 0] - E[n-2, 0])/3
        E[n,Nr-1] = (2j*dZ*(dEdr2L(E, n) + R[-1]*(E[n-1, -2] - E[n-1, -1]/dR) + E3(E, -1, n)) + 4*E[n-1, -1] - E[n-2, -1])/3

        if n % 100 == 0:
            print(n)
            print(trapz(np.abs(E[n, :])))



    np.save("hest", E)
    plot(R, np.abs(E[0, :]))
    # plot(R, E[1, :])
    # plot(R, E[200, :])
    # plot(R, E[300, :])
    # plot(R, E[400, :])
    plot(R, np.abs(E[40000, :]))
    plot(R, np.abs(E[-1, :]))




    print(np.trapz(np.abs(E[0, :])))
    print(np.trapz(np.abs(E[-1, :])))
    # plot(R, E[-1, :])
    
    # if (plot3d):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')
    #
    #
    #     X,Y = np.meshgrid(R,np.linspace(0,Zmax,Nz+1));
    #
    #     ax.plot_surface(X,Y,np.abs(E)**2);
    
    
    show()



solve_PDE2()