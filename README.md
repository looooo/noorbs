# noorbs

This is a little library (c++ and python-bindings) to compute 1D, 2D, 3D non uniform rational bsplines. Computing derivatives analytically can be used to solve pde's.

## Dependencies  
To build the library cmake, eigen, pybind11 are necessary. Numpy is essential to use the library from python.

## Usage in python:  

```python
import numpy as np
import nurbs

# setup the basic input for the basis
u_knots = nurbs.create_knots_vector(u_min=0, u_max=1, degree=2, num_poles=4)
v_knots = nurbs.create_knots_vector(u_min=0, u_max=1, degree=3, num_poles=4)
u_degreee = 1
v_degree = 1
weights = np.array([1.] * 8)

# create the nurbs-object
a = nurbs.NurbsBase2D(u_knots, v_knots, weights, u_degreee, v_degree)

# precompute the derivatives (setup functions)
a.computeFirstDerivatives()

# specify points where we want to evaluate the nurbs-surface
u = np.linspace(u_knots[0], u_knots[-1], 21)
v = np.linspace(v_knots[0], v_knots[-1], 21)
uv = np.array([[i, j] for i in u for j in v])

# compute influence-matrix and derivatives at given points
mat = a.getInfluenceMatrix(uv)
dv = a.getDVMatrix(uv)
du = a.getDUMatrix(uv)

# now we can easily compute the interpolation-data by matrix-products
data = uv[0] ** 2 + uv[1] ** 2 
int_data = mat @ data
int_data_du = du @ data
int_data_dv = dv @ data
```


Still there are many things missing like:  
- higher derivatives (theoretically possible)
- faster computation of influence and derivatives
- ...