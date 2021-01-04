import gmsh
from pyelmer import gmsh_utils

factory = gmsh.model.occ

class Parameters:
    pass


class Model:
    def __init__(self):
        self.shapes_3d = []
        self.shapes_2d = []
        self.shapes_1d = []
        self.shapes_0d = []

    def __enter__(self):
        gmsh.initialize()
        gmsh.model.add("occ model")
        return self
    
    def __exit__(self, *args):
        gmsh.finalize()

    def show(self):
        gmsh.fltk.run()

    def add_shape(self, shape, name):
        if shape.dim == 3:
            self.shapes_3d.append(shape)
        elif shape.dim == 2:
            self.shapes_2d.append(shape)
        elif shape.dim == 1:
            self.shapes_1d.append(shape)
        elif shape.dim == 0:
            self.shapes_0d.append(shape)
    
    def make_physical(self, dimension):
        if dimension == 3:
            for shape in self.shapes_3d:
                shape.make_physical()
        if dimension == 2:
            for shape in self.shapes_2d:
                shape.make_physical()
        if dimension == 1:
            for shape in self.shapes_1d:
                shape.make_physical()
        if dimension == 0:
            for shape in self.shapes_0d:
                shape.make_physical()

    def synchronize(self):
        factory.synchronize()


class Shape:
    def __init__(self, model, dim, name, geo_ids=[]):
        self.dim = dim
        self.name = name
        self.geo_ids = geo_ids
        self.params = Parameters()
        self.ph_id = -1
        self.mesh_size = 0
        self.neighbors = []
        model.add_shape(self, self.name)

    @property
    def dimtags(self):
        return [(self.dim, x) for x in self.geo_ids]

    @property
    def boundaries(self):
        bndry = []
        for geo_id in self.geo_ids:
            bndry += gmsh_utils.get_boundaries(self.dim, geo_id)
        return [x for x in bndry if bndry.count(x) == 1]  # exclude boundaries inside shape

    def set_interface(self, shape):
        self.neighbors.append(shape)
        factory.fragment(self.dimtags, shape.dimtags)
        factory.synchronize()

    def get_interface(self, shape):
        own = self.boundaries
        other = shape.boundaries
        return [x for x in own if x in other]

    def make_physical(self):
        self.ph_id = gmsh_utils.add_physical_group(self.dim, self.geo_ids, self.name)


# class Boundary:
#     def __init__(self, dim, geo_ids, name):
#         self.dim = dim
#         self.geo_ids = geo_ids
#         self.name = name
#         self.ph_id = -1
    
#     def make_physical(self):
#         self.ph_id = gmsh_utils.add_physical_group(self.dim, self.geo_ids, self.name)

