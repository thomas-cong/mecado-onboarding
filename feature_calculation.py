import json
import vtk
from OCC.Core.BRepTools import breptools
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.BRep import BRep_Tool

from VTKTool import BRepMesh




def mean_curvature(brep: BRepMesh):
    '''
    Calculate mean curvature for each face
    '''
    result = {}
    for face_id, _, face_points in brep.iter_faces():
        face = brep.face_shapes.get(face_id)
        if face is None:
            continue
        # Parametric bounds on the face
        umin, umax, vmin, vmax = breptools.UVBounds(face)
        u_mid = 0.5 * (umin + umax)
        v_mid = 0.5 * (vmin + vmax)

        # Underlying surface and curvature props
        surf = BRep_Tool.Surface(face)
        props = GeomLProp_SLProps(surf, u_mid, v_mid, 2, 1.0e-6)

        if props.IsCurvatureDefined():
            kmax = props.MaxCurvature()
            kmin = props.MinCurvature()
            H = 0.5 * (kmax + kmin)  # mean curvature
        else:
            H = 0.0  # planar or undefined curvature

        result[face_id] = H
        
    return result

if __name__ == "__main__":
    Brep = BRepMesh("./rearrimbolt.stp", 0.1, 0.5)
    print(mean_curvature(Brep))