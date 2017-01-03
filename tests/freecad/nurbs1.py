import numpy as np
import nurbs

u_knots = np.array([ 0, 0, 0, 1, 1, 1], float)
v_knots = np.array([ 0, 0, 0, 1, 1, 1], float)
u_degreee = 2
v_degree = 2
weights = np.array([1., 1., 1., 1., 1, 1, 1, 1, 1], float)

u_dess = 20
v_dess = 20

a = nurbs.NurbsBase2D(u_knots, v_knots, weights, u_degreee, v_degree)
a.computeFirstDerivatives()
u = np.linspace(u_knots[0], u_knots[-1], u_dess)
v = np.linspace(v_knots[0], v_knots[-1], v_dess)

uv = [[i, j] for i in u for j in v]

x_points = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
y_points = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
z_points = np.array([0, 0, 0, 0, 10, 0, 0, 0, 0])
pos = a.getInfluenceMatrix(np.array(uv)) * np.array([x_points, y_points, z_points]).T

a.getDuVector(np.array([0, 0]))
du = a.getDuMatrix(np.array(uv)) * z_points
dv = a.getDvMatrix(np.array(uv)) * z_points
du = du ** 2 + dv ** 2

colors = (du - min(du))/ (max(du) - min(du))
colors = np.array([colors, np.zeros(len(colors)), np.zeros(len(colors))]).T.tolist()
from pivy import coin
# # from pivy_primitives_new_new import Marker

result = coin.SoSeparator()

myBinding = coin.SoMaterialBinding()
myBinding.value = coin.SoMaterialBinding.PER_VERTEX_INDEXED
result.addChild(myBinding)

myMaterial = coin.SoMaterial()
myMaterial.diffuseColor.setValues(0, len(colors), colors)
result += myMaterial

# Define coordinates for vertices
myCoords = coin.SoCoordinate3()
myCoords.point.setValues(0, len(pos), pos.tolist())
result += myCoords

# Define the QuadMesh.
myQuadMesh = coin.SoQuadMesh()
myQuadMesh.verticesPerRow = u_dess
myQuadMesh.verticesPerColumn = v_dess

result += myQuadMesh

import FreeCADGui
sg = FreeCADGui.ActiveDocument.ActiveView.getSceneGraph()
sg += result

