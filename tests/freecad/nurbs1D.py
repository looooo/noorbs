import nurbs
import numpy as np
from matplotlib import pyplot as plt

points = np.array([0, 1, 0], float)
weights = np.array([1, 1, 1], float)
knots =  np.array([0, 0, 0, 1, 1, 1], float)
t_sequence = np.linspace(0, 1, 100)
degree = 2

base = nurbs.NurbsBase1D(knots, weights, degree)
base.computeFirstDerivatives()
print(base.getDuVector(0.4))

plt.plot(base.getInfluenceMatrix(t_sequence).toarray())
plt.show()

plt.plot(t_sequence, base.getInfluenceMatrix(t_sequence) * points)
plt.show()

plt.plot(base.getDuMatrix(t_sequence) * points)
plt.show()