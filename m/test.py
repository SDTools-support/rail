
import gmsh

occ = gmsh.model.occ
try:
    gmsh.initialize()
    gmsh.open("Rail_UIC60.stp")
    print(gmsh.model.getEntities())
    occ.removeAllDuplicates()
    occ.synchronize()   
    occ.remove([(1, 1)])
    occ.synchronize()    
    print("all entities\n", gmsh.model.getEntities())
    new_line = occ.addLine(37, 2, tag=1)
    occ.synchronize()
    
    lines = [e[1] for e in gmsh.model.getEntities(dim=1)]
    print("lines\n", lines)
    cloop = occ.addCurveLoop(lines)
    surf = occ.addPlaneSurface([cloop])
    occ.extrude([(2,surf)], 20, 0, 0)
    occ.synchronize()
    gmsh.model.mesh.setSize([(0, 8),(0, 7),(0, 6),(0, 43),(0, 44),(0, 45)], 1)
    occ.synchronize()
    gmsh.model.mesh.setOrder(1)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 5)
    gmsh.model.mesh.generate(3)
    gmsh.fltk.run()
    
finally:
    gmsh.finalize()