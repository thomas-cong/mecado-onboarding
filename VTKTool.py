from vtkmodules.vtkCommonCore import vtkPoints, vtkIntArray, vtkStringArray
from vtkmodules.vtkCommonDataModel import vtkCellArray
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopLoc import TopLoc_Location
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
import json
import vtkmodules.vtkRenderingFreeType  # side-effect import also works
import vtkmodules.vtkRenderingOpenGL2
import vtk

import feature_calculation

# Visualizer Imports
from vtkmodules.vtkRenderingCore import (
    vtkPolyDataMapper, vtkActor, vtkRenderer, vtkRenderWindow,
    vtkRenderWindowInteractor, vtkTextActor
)
from vtkmodules.vtkRenderingCore import vtkTextProperty
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkRenderingCore import vtkCellPicker
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

class BRepMesh:
    def __init__(self, step_path: str, linear_deflection: float = 0.1, angular_deflection: float = 0.5, quant_scale: int = 1_000_000):
        # load and mesh
        self.shape = self._load_step(step_path)
        self.mesh = self._mesh_shape(self.shape, linear_deflection, angular_deflection)

        # storage
        self.pts = vtkPoints()
        self.cells = vtkCellArray()
        self.face_ids = vtkIntArray()
        self.face_ids.SetName("FaceID")
        self.face_labels = vtkStringArray()
        self.face_labels.SetName("FaceLabel")
        self.point_map = {}  # quantized coord -> global point id
        self.quant_scale = quant_scale
        self.step_path = step_path
        self.face_shapes = {}
        self.build()

    def _load_step(self, path: str):
        '''
        Load the STEP file
        '''
        reader = STEPControl_Reader()
        reader.ReadFile(path)
        reader.TransferRoots()
        return reader.OneShape()

    def _mesh_shape(self, shape, linear_deflection: float, angular_deflection: float):
        '''
        Mesh the shape
        '''
        mesh = BRepMesh_IncrementalMesh(shape, linear_deflection, False, angular_deflection, True)
        mesh.Perform()
        return mesh

    def _quantize(self, p):
        '''
        Quantize a point (3D helper)
        '''
        s = self.quant_scale
        return (int(round(p.X() * s)), int(round(p.Y() * s)), int(round(p.Z() * s)))

    def build(self):
        '''
        Build the mesh and fill arrays, labels, etc.
        '''
        self._enumerate_faces(self._add_face_to_arrays)
        self.built = True
        try:
            json_path = self.step_path.replace(".stp", ".json")
            self.label_faces_from_json(json_path)
        except:
            print("failed to find json")
            pass

    def _enumerate_faces(self, callback):
        '''
        Enumerate all faces in the mesh
        '''
        exp = TopExp_Explorer(self.shape, TopAbs_FACE)
        face_id = 0
        while exp.More():
            face_id += 1
            face = topods.Face(exp.Current())
            self.face_shapes[face_id] = face
            loc = TopLoc_Location()
            tri = BRep_Tool.Triangulation(face, loc)
            callback(face, face_id, tri, loc)
            exp.Next()

    def _add_face_to_arrays(self, face, face_id, tri, loc):
        '''
        Add a face to the arrays
        '''
        if tri is None:
            return

        T = loc.Transformation()

        # face-local node id -> global point id
        local_to_global = {}
        for i in range(1, tri.NbNodes() + 1):
            gp = tri.Node(i)
            gp.Transform(T)
            key = self._quantize(gp)

            if key in self.point_map:
                gid = self.point_map[key]
            else:
                gid = self.pts.GetNumberOfPoints()
                self.pts.InsertNextPoint(gp.X(), gp.Y(), gp.Z())
                self.point_map[key] = gid

            local_to_global[i] = gid

        # triangles using global ids
        for t_i in range(1, tri.NbTriangles() + 1):
            n1, n2, n3 = tri.Triangle(t_i).Get()
            self.cells.InsertNextCell(3)
            self.cells.InsertCellPoint(local_to_global[n1])
            self.cells.InsertCellPoint(local_to_global[n2])
            self.cells.InsertCellPoint(local_to_global[n3])
            self.face_ids.InsertNextValue(face_id)
    def label_faces_from_json(self, path: str):
        '''
        Label faces from a json file
        '''
        self.check_built()
        with open(path, "r") as f:
            data = json.load(f)
        id_to_label = {int(item["name"]): item["label"] for item in data["Faces"]}
        n_cells = self.face_ids.GetNumberOfTuples()
        for i in range(n_cells):
            face_id = int(self.face_ids.GetValue(i))
            label = id_to_label.get(face_id, "")
            self.face_labels.InsertNextValue(label)
    def make_polydata(self):
        '''
        Make a vtkPolyData object from the mesh
        '''
        self.check_built()
        poly = vtkPolyData()
        poly.SetPoints(self.pts)
        poly.SetPolys(self.cells)
        # face_ids is per-cell data
        poly.GetCellData().AddArray(self.face_ids)
        poly.GetCellData().AddArray(self.face_labels)
        return poly

    def write_vtp(self, path: str):
        '''
        Write the mesh to a VTP file
        '''
        self.check_built()
        writer = vtkXMLPolyDataWriter()
        writer.SetFileName(path)
        writer.SetInputData(self.make_polydata())
        writer.Write()
    def check_built(self):
        '''
        Check that the mesh has been built and lists are filled
        '''
        if not self.built:
            raise ValueError("BRepMesh must be built before this operation")
    def iter_faces(self):
        '''
        Iterate over all faces in the mesh and return triangles
        '''
        self.check_built() # check that lists are filled 
        id_list = vtk.vtkIdList() # creates empty object of vtkIdList
        self.cells.InitTraversal() # resets traversal to start from beginning
        n_cells = self.face_ids.GetNumberOfTuples() # determines how many unique faces there are (which is also how many cells)
        for i in range(n_cells): # loop over all cells
            ok = self.cells.GetNextCell(id_list) # gets next cell and stores it in id_list
            if not ok:
                raise ValueError("Failed to get cell") # if can't get cell fail
            num_ids = id_list.GetNumberOfIds() # get number of ids in this cell (three)
            points = [self.pts.GetPoint(id_list.GetId(j)) for j in range(num_ids)] # get points for this cell
            face_id = self.face_ids.GetValue(i) # get face id
            face_label = self.face_labels.GetValue(i) # get face label
            yield (face_id, face_label, points) # yield face id, label, and points
    def calculate_features(self):
        '''
        Calculate features for each face
        '''
        self.check_built()
        features = {}
        curvatures = feature_calculation.curvatures(self)
        normals = feature_calculation.normals(self)
        for face_id, face_label, points in self.iter_faces():
            features[face_id] = {
                "label": face_label,
                "features": {
                    "curvatures": curvatures[face_id],
                    "normals": normals[face_id]
                }
            }
        return features
    def write_json(self):
        '''
        Write the features to a json file
        '''
        self.check_built()
        features = self.calculate_features()
        formatted_features = {"Faces": []}
        for face_id, feature_values in features.items():
            formatted_features["Faces"].append({"name": face_id, **feature_values})
        with open(self.step_path.replace(".stp", ".json"), "w") as f:
            json.dump(formatted_features, f)

class VTPVisualizer:
    def __init__(self, path: str, brep: BRepMesh | None = None):
        self.reader = vtkXMLPolyDataReader()
        self.reader.SetFileName(path)
        self.reader.Update()
        self.polydata = self.reader.GetOutput()
        self.brep = brep
        cd = self.polydata.GetCellData()
        print("CellData arrays:", cd.GetNumberOfArrays())
        for i in range(cd.GetNumberOfArrays()):
            print("CellData array", i, ":", cd.GetArrayName(i))
        # Create a per-cell color array: fillet -> green, others -> gray
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("FaceColor")
        n_cells = self.polydata.GetNumberOfCells()
        arr_labels = cd.GetAbstractArray("FaceLabel")
        for i in range(n_cells):
            label = arr_labels.GetValue(i) if arr_labels is not None else ""
            if label == "fillet":
                colors.InsertNextTuple3(0, 255, 0)   # green
            else:
                colors.InsertNextTuple3(200, 200, 200)  # gray
        cd.SetScalars(colors)

        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputData(self.polydata)
        self.mapper.SetScalarModeToUseCellData()
        self.mapper.SetColorModeToDirectScalars()

        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        
        self.renderer = vtkRenderer()
        self.renderer.AddActor(self.actor)

        self.render_window = vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)

        self.render_window_interactor = vtkRenderWindowInteractor()
        self.render_window_interactor.SetRenderWindow(self.render_window)

        style = vtkInteractorStyleTrackballCamera()
        self.render_window_interactor.SetInteractorStyle(style)

    def enable_hover_face_labels(self):
        self.text_actor = vtkTextActor()

        prop = self.text_actor.GetTextProperty()
        prop.SetFontSize(18)
        prop.SetColor(0, 0, 0)
        prop.SetBackgroundColor(1.0, 1.0, 1.0)
        prop.SetBackgroundOpacity(1.0)
        self.text_actor.SetDisplayPosition(10, 10)
        self.text_actor.SetInput("FaceID: -")
        self.renderer.AddActor2D(self.text_actor)

        self.picker = vtkCellPicker()
        self.render_window_interactor.SetPicker(self.picker)

        ren = self.renderer
        ren_win = self.render_window
        ren_win_interactor = self.render_window_interactor
        poly = self.polydata
        text_actor = self.text_actor
        picker = self.picker

        def on_mouse_move(obj, event):
            x,y = ren_win_interactor.GetEventPosition()
            if not picker.Pick(x, y, 0, ren):
                text_actor.SetInput("FaceID: -")
                ren_win.Render()
                return
            cell_id = picker.GetCellId()
            if cell_id < 0:
                text_actor.SetInput("FaceID: -")
                ren_win.Render()
                return
            arr_ids = poly.GetCellData().GetArray("FaceID")
            arr_labels = poly.GetCellData().GetAbstractArray("FaceLabel")
            if arr_ids is None or arr_labels is None:
                text_actor.SetInput("FaceID: (none)")
                ren_win.Render()
                return
            face_id = int(arr_ids.GetValue(cell_id))
            face_label = arr_labels.GetValue(cell_id)   
            text_actor.SetInput(f"FaceID: {face_id}\nFaceLabel: {face_label}")
            ren_win.Render()
        self.render_window_interactor.AddObserver("MouseMoveEvent", on_mouse_move)

        def on_left_button_press(obj, event):
            x, y = ren_win_interactor.GetEventPosition()
            if not picker.Pick(x, y, 0, ren):
                return
            cell_id = picker.GetCellId()
            if cell_id < 0:
                return

            arr_ids = poly.GetCellData().GetArray("FaceID")
            arr_labels = poly.GetCellData().GetAbstractArray("FaceLabel")
            if arr_ids is None or arr_labels is None:
                return

            face_id = int(arr_ids.GetValue(cell_id))
            current_label = arr_labels.GetValue(cell_id)

            # Toggle between no label (non-fillet) and 'fillet'
            new_label = "fillet" if current_label == "" else ""

            # Update all cells in the polydata that belong to this face_id
            n_cells = poly.GetNumberOfCells()
            colors = poly.GetCellData().GetScalars()
            for cid in range(n_cells):
                if int(arr_ids.GetValue(cid)) == face_id:
                    arr_labels.SetValue(cid, new_label)
                    if new_label == "fillet":
                        colors.SetTuple3(cid, 0, 255, 0)
                    else:
                        colors.SetTuple3(cid, 200, 200, 200)

            # Keep the underlying BRepMesh labels in sync and rewrite JSON/features
            if self.brep is not None:
                brep_ids = self.brep.face_ids
                brep_labels = self.brep.face_labels
                n_brep_cells = brep_ids.GetNumberOfTuples()
                for i in range(n_brep_cells):
                    if int(brep_ids.GetValue(i)) == face_id:
                        brep_labels.SetValue(i, new_label)
                # Recompute features and write updated JSON
                self.brep.write_json()

            # Notify VTK that scalars and polydata have changed so colors update
            if colors is not None:
                colors.Modified()
                poly.GetCellData().SetScalars(colors)
            poly.Modified()

            text_actor.SetInput(f"FaceID: {face_id}\nFaceLabel: {new_label}")
            ren_win.Render()

        self.render_window_interactor.AddObserver("LeftButtonPressEvent", on_left_button_press)
    def start(self):
        self.renderer.ResetCamera()
        self.render_window.Render()
        self.render_window_interactor.Initialize()
        self.render_window_interactor.Start()
if __name__ == "__main__":
    brep = BRepMesh("./rearrimbolt.stp", 0.1, 0.5)
    brep.write_json()
    brep.write_vtp("./rearrimbolt.vtp")
    viz = VTPVisualizer("./rearrimbolt.vtp", brep=brep)
    viz.enable_hover_face_labels()
    viz.start()