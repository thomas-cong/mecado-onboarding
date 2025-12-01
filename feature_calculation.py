import json
import vtk
from OCC.Core.BRepTools import breptools
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.BRep import BRep_Tool
import numpy as np



def curvatures(brep):
    '''
    Calculate min, max, mean curvature for each face
    flat: 0 min and max curvature
    fillet: non-zero min and max curvature
    '''
    result = {}
    for face_id, _, face_points in brep.iter_faces():
        face = brep.face_shapes.get(face_id)
        if face is None:
            continue
        # Parametric bounds on the face
        umin, umax, vmin, vmax = breptools.UVBounds(face)

        # Underlying surface
        surf = BRep_Tool.Surface(face)

        # Sample a small UV grid over the face and average curvatures
        n_u, n_v = 3, 3
        du = (umax - umin) / (n_u - 1) if n_u > 1 else 0.0
        dv = (vmax - vmin) / (n_v - 1) if n_v > 1 else 0.0

        sum_kmin = 0.0
        sum_kmax = 0.0
        sum_H = 0.0
        n_valid = 0

        for iu in range(n_u):
            u = umin + iu * du
            for iv in range(n_v):
                v = vmin + iv * dv
                props = GeomLProp_SLProps(surf, u, v, 2, 1.0e-6)
                if not props.IsCurvatureDefined():
                    continue
                kmax = props.MaxCurvature()
                kmin = props.MinCurvature()
                H = 0.5 * (kmax + kmin)
                sum_kmin += kmin
                sum_kmax += kmax
                sum_H += H
                n_valid += 1

        if n_valid > 0:
            kmin_avg = sum_kmin / n_valid
            kmax_avg = sum_kmax / n_valid
            H_avg = sum_H / n_valid
        else:
            # treat as planar / undefined curvature
            kmin_avg = 0.0
            kmax_avg = 0.0
            H_avg = 0.0

        result[face_id] = (kmin_avg, kmax_avg, H_avg)
    return result
def normals(brep):
    '''
    Calculate max, min, variance, mean of dot products of normals for each face
    flat: 0 variance
    fillet: non-zero variance
    '''
    normals = {}
    def calculate_normals(face_points):
        if len(face_points) < 3:
            return None
        v1 = np.array(face_points[1]) - np.array(face_points[0])
        v2 = np.array(face_points[2]) - np.array(face_points[0])
        n = np.cross(v1, v2)
        n = n / np.linalg.norm(n)
        return n
    for face_id, _, face_points in brep.iter_faces():
        n = calculate_normals(face_points)
        if n is None:
            continue
        normals[face_id] = normals.get(face_id, [])
        normals[face_id].append(n)
    
    result = {}
    for face_id, n_list in normals.items():
        mean_n = np.mean(n_list, axis=0) # unit norm
        mean_norm = np.linalg.norm(mean_n)
        if mean_norm == 0.0:
            result[face_id] = 1.0  # or some max variance
            continue
        mean_n = mean_n / mean_norm
        dots = np.array([np.dot(n, mean_n) for n in n_list])
        max_dot = np.max(dots)
        min_dot = np.min(dots)
        var = np.var(dots)
        mean_dot = np.mean(dots)
        result[face_id] = (max_dot, min_dot, var, mean_dot) # variance of dot products
    return result
    

if __name__ == "__main__":
    from VTKTool import BRepMesh

    Brep = BRepMesh("./rearrimbolt.stp", 0.1, 0.5)
    print(curvatures(Brep))
    # print(normals(Brep))