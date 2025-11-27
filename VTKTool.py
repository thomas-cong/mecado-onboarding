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

    def _load_step(self, path: str):
        reader = STEPControl_Reader()
        reader.ReadFile(path)
        reader.TransferRoots()
        return reader.OneShape()

    def _mesh_shape(self, shape, linear_deflection: float, angular_deflection: float):
        mesh = BRepMesh_IncrementalMesh(shape, linear_deflection, False, angular_deflection, True)
        mesh.Perform()
        return mesh

    def _quantize(self, p):
        s = self.quant_scale
        return (int(round(p.X() * s)), int(round(p.Y() * s)), int(round(p.Z() * s)))

    def build(self):
        self._enumerate_faces(self._add_face_to_arrays)
        try:
            json_path = self.step_path.replace(".stp", ".json")
            self.label_faces_from_json(json_path)
        except:
            print("failed to find json")
            pass

    def _enumerate_faces(self, callback):
        exp = TopExp_Explorer(self.shape, TopAbs_FACE)
        face_id = 0
        while exp.More():
            face_id += 1
            face = topods.Face(exp.Current())
            loc = TopLoc_Location()
            tri = BRep_Tool.Triangulation(face, loc)
            callback(face, face_id, tri, loc)
            exp.Next()

    def _add_face_to_arrays(self, face, face_id, tri, loc):
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
        with open(path, "r") as f:
            data = json.load(f)
        id_to_label = {int(item["name"]): item["label"] for item in data["Faces"]}
        n_cells = self.face_ids.GetNumberOfTuples()
        for i in range(n_cells):
            face_id = int(self.face_ids.GetValue(i))
            label = id_to_label.get(face_id, "")
            self.face_labels.InsertNextValue(label)
    def make_polydata(self):
        poly = vtkPolyData()
        poly.SetPoints(self.pts)
        poly.SetPolys(self.cells)
        # face_ids is per-cell data
        poly.GetCellData().AddArray(self.face_ids)
        poly.GetCellData().AddArray(self.face_labels)
        return poly

    def write_vtp(self, path: str):
        writer = vtkXMLPolyDataWriter()
        writer.SetFileName(path)
        writer.SetInputData(self.make_polydata())
        writer.Write()
class VTPVisualizer:
    def __init__(self, path: str):
        self.reader = vtkXMLPolyDataReader()
        self.reader.SetFileName(path)
        self.reader.Update()
        self.polydata = self.reader.GetOutput()
        cd = self.polydata.GetCellData()
        print("CellData arrays:", cd.GetNumberOfArrays())
        for i in range(cd.GetNumberOfArrays()):
            print("CellData array", i, ":", cd.GetArrayName(i))
        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputData(self.polydata)

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
    def start(self):
        self.renderer.ResetCamera()
        self.render_window.Render()
        self.render_window_interactor.Initialize()
        self.render_window_interactor.Start()
if __name__ == "__main__":
    tool = BRepMesh("./rearrimbolt.stp", 0.1, 0.5)
    tool.build()
    tool.write_vtp("./rearrimbolt.vtp")
    viz = VTPVisualizer("./rearrimbolt.vtp")
    viz.enable_hover_face_labels()
    viz.start()