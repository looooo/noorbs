import nurbs
import Mesh

obj = Gui.Selection.getSelection()[0] # obj must be a Mesh (Mesh-Design->Meshes->Create-Mesh)
mesh = Mesh.Mesh(obj.Mesh) # copy of the mesh to set new vertices later on
points = [[i.x, i.y, i.z] for i in obj.Mesh.Points]
faces = [list(i) for i in  obj.Mesh.Topology[1]]
flatten = nurbs.LscmRelax(points, faces, [])
flatten.lscm()     # least square conformal map
flatten.relax(0.9) # fem step value 0-0.999 (for double curved surfaces call this function x times)
for i, point in enumerate(flatten.flat_vertices):
    mesh.setPoint(i, App.Vector(*point))
Mesh.show(mesh)
