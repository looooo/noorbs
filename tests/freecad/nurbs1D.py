import nurbs
import numpy as np
from matplotlib import pyplot as plt

points = np.array([0, 1, 0, 1, 0, 1], float)
weights = np.array([1, 1, 1, 1, 1, 1], float)
knots =  np.array([0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1], float)
t_sequence = np.linspace(0, 1, 1000)
degree = 2

base = nurbs.NurbsBase1D(knots, weights, degree)
# print(base.getInfluenceMatrix(t_sequence) * points)

plt.plot(base.getInfluenceMatrix(t_sequence).toarray())
plt.show()


plt.plot(t_sequence, base.getInfluenceMatrix(t_sequence) * points)
plt.show()