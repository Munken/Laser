from pylab import *

X = np.linspace(-np.pi, np.pi, 256,endpoint=True)
C,S,G = np.cos(X), np.sin(X), np.exp(-(np.square(X)))

plot(X,C)
plot(X,S)
plot(X,G)

show()