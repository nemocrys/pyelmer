import gmsh
import numpy as np
from pyelmer import gmsh_utils

factory = gmsh.model.occ
field = gmsh.model.mesh.field

#TODO point select of geometries / boundaries / ...?

class GmshError(Exception):
    pass


class GeometryError(Exception):
    pass


class Parameters:
    """Dummy class, used as container to store geometry parameters of
    the shapes.
    """
    def __init__(self):
        self.T_init = 0


class Model:
    """Wrapper for Gmsh model. This class is used to manage the shapes
    and physical groups, and provides high-level access to some major
    functionalities of the gmsh API.
    """
    def __init__(self, name='model'):
        """Create gmsh model.

        Args:
            name (str, optional): Name of the model. Defaults to 'model'.
        """
        self._shapes = []
        self.mesh_restrictions = []
        self.min_field = -1
        self._physical = False

        gmsh.initialize()
        gmsh.model.add(name)  
    
    def __del__(self):
        gmsh.finalize()

    def __getitem__(self, name):
        """Get shape that is part of this model.

        Args:
            name (str): Name of the shape

        Returns:
            Shape: shape object with shape.name == name.
        """
        for shape in self._shapes:
            if shape.name == name:
                return shape
        raise GeometryError(f'Shape {name} does not exist.')
    
    def __repr__(self):
        shapes = [s.name for s in self._shapes]
        return f'Gmsh model created with pyelmer.\nShapes: {shapes}'

    def show(self):
        """Run gmsh GUI."""
        gmsh.fltk.run()

    def _add_shape(self, shape):
        """Add a shape to the model.

        Args:
            shape (Shape): shape object containing the information about
                           the geometry
        """
        self._shapes.append(shape)

    def _apply_restrictions(self):
        self.min_field = field.add('Min')
        field.setNumbers(self.min_field, 'FieldsList', [x.field for x in self.mesh_restrictions])
        field.setAsBackgroundMesh(self.min_field)

    def deactivate_characteristic_length(self):
        gmsh.option.setNumber('Mesh.CharacteristicLengthFromPoints', 0)
        gmsh.option.setNumber('Mesh.CharacteristicLengthFromCurvature', 0)
        gmsh.option.setNumber('Mesh.CharacteristicLengthExtendFromBoundary', 0)

    def set_characteristic_length(self, char_length, dimensions=[0]):
        for dim in dimensions:
            gmsh.model.mesh.setSize(gmsh.model.getEntities(dim), char_length)

    def generate_mesh(self, dimension=2, order=1, size_factor=1, smoothing=1):
        self._apply_restrictions()
        gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', size_factor)
        gmsh.option.setNumber('Mesh.Smoothing', smoothing)
        gmsh.model.mesh.generate(dimension)
        gmsh.model.mesh.setOrder(order)    

    def get_shapes(self, dimension, name=''):
        """Get shapes of certain dimension, optionally filtered by name.

        Args:
            dimension (int): Dimension of shapes
            name (str, optional): Only shapes with this string in their
                                  name.

        Returns:
            list: shape objects with shape.dim = dimension and name in
                  shape.name
        """
        return [s for s in self._shapes if s.dim == dimension and name in s.name]

    def make_physical(self):
        """Convert all shapes into physical groups.

        Raises:
            GmshError: If this function is called more than once.
        """
        if self._physical:
            raise GmshError('This model is already physical.')
        for shape in self._shapes:
            shape._make_physical()

    def synchronize(self):
        """Synchronize gmsh geometry kernel.
        """
        # if self._physical:
        #     raise GmshError('The model is physical. Synchronizing would break it.')
        factory.synchronize()

    def remove_shape(self, shape):
        """Remove shape from model.

        Args:
            shape (Shape): shape object to be removed
        """
        self._shapes.remove(shape)

    def write_msh(self, file_name):
        gmsh.write(file_name)

    def  set_const_mesh_sizes(self):
        for shape in self._shapes:
            if shape.mesh_size == 0:
                print(f'Warning: Mesh size = 0 for {shape.name}. Ignoring this shape...')
            else:
                print(shape.name, shape.mesh_size)
                MeshControlConstant(self, shape.mesh_size, shapes=[shape])

class Shape:
    """Wrapper for any kind of shape, that shall be part of the final
    model. Shapes may be 2 or 3D objects, lines or points.
    A shape may consist of multiple gmsh entities of the same dimension.
    """
    def __init__(self, model, dim, name, geo_ids=[]):
        """Create a shape.

        Args:
            model (Model): Gmsh model to which the shape will be added
            dim (int): Dimension of the shape
            name (str): Name of the shape
            geo_ids (list, optional): Geometry IDs of the shape.
                                      They may be added later.
        """
        self.model = model
        self.dim = dim
        self.name = name
        self.geo_ids = geo_ids
        self.params = Parameters()
        self.ph_id = -1
        self.mesh_size = 0
        model._add_shape(self)

    def __iadd__(self, other):
        if type(other) is list:
            self.geo_ids += [x for x in other if x not in self.geo_ids]
        if type(other) is int:
            self.geo_ids.append(other)
        if type(other) is Shape:
            self.geo_ids += [x for x in other.geo_ids if x not in self.geo_ids]
            self.model.remove_shape(other)
        return self

    def __isub__(self, other):
        if type(other) is list:
            self.geo_ids = [x for x in self.geo_ids if x not in other]
        if type(other) is int and other in self.geo_ids:
            self.geo_ids.remove(other)
        if type(other) is Shape:
            self.geo_ids = [x for x in self.geo_ids if x not in other.geo_ids]
        return self

    @property
    def dimtags(self):
        """Gmsh dimension tags.

        Returns:
            list: Gmsh dim-tags of entities in shape.
        """
        return [(self.dim, x) for x in self.geo_ids]

    @property
    def boundaries(self):
        """Boundaries of the shape (dimension: dimension of shape - 1).

        Returns:
            list: Tags of the external boundaries of the shape.

        Note: Internal boundaries, e.g. a point in between two lines
        that form a shape, are excluded. Only the endpoints of the lines
        are cosidered as boundaries.
        """
        bndry = []
        for geo_id in self.geo_ids:
            bndry += gmsh_utils.get_boundaries(self.dim, geo_id)
        return [x for x in bndry if bndry.count(x) == 1]  # exclude boundaries inside shape

    @property
    def bounding_box(self):
        """Get the bounding box of this shape.

        Returns:
            list[float]: [x_min, y_min, z_min, x_max, y_max, z_max]
        """
        boxes = np.array([factory.getBoundingBox(self.dim, tag) for tag in self.geo_ids])
        return [boxes[:, 0].min(), boxes[:, 1].min(), boxes[:, 2].min(),
                boxes[:, 3].max(), boxes[:, 4].max(), boxes[:, 4].max()]

    def set_interface(self, shape):
        """Get to know the other shape and remove duplicate boundaries.
        Only required if shapes are in contact with each other.

        Args:
            shape (Shape): Other shape that is in contact with this.
        """
        # self.neighbors.append(shape)
        # shape.neighbors.append(self)
        factory.fragment(self.dimtags, shape.dimtags)
        factory.synchronize()

    def get_interface(self, shape):
        """Get boundaries that lay in between two shapes.

        Args:
            shape (Shape): Other shape.

        Returns:
            list: Tags of boundary elements.
        """
        own = self.boundaries
        other = shape.boundaries
        return [x for x in own if x in other]

    def get_boundaries_in_box(self, x, y, z=[0, 0], eps=1e-6, one_only=False):
        """Get boundaries of the shape with a box-select. Only
        boundaries that are completely inside the box are returned.

        Args:
            x (list): [x_min, x_max] limits for x-coordinate
            y (list): [y_min, y_max] limits for y-coordinate
            z (list, optional): [z_min, z_max] limits for z coordinate.
            eps (float, optional): Sensitivity. Defaults to 1e-6.
            one_only (bool, optional): Search for only one boundary;
            raise GeometryError if more / less boundaries are found.

        Returns:
            list: Boundary tags.
            Integer return if only_one was set to true.
        """
        dimtags = gmsh.model.getEntitiesInBoundingBox(x[0] - eps, y[0] - eps, z[0] - eps,
                                                      x[1] + eps, y[1] + eps, z[1] + eps,
                                                      self.dim - 1)
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.boundaries]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(f'Found {len(tags_filtered)} instead of only one boundary.')
            else:
                return tags_filtered[0]
        else:
            return tags_filtered

    def get_part_in_box(self, x, y, z=[0, 0], eps=1e-6, one_only=False):
        """Get part of the shape within a box. Only the parts completely
        inside the box are returned

        Args:
            x (list): [x_min, x_max] limits for x-coordinate
            y (list): [y_min, y_max] limits for y-coordinate
            z (list, optional): [z_min, z_max] limits for z coordinate.
            eps (float, optional): Sensitivity. Defaults to 1e-6.
            one_only (bool, optional): Search for only one part;
            raise GeometryError if more / less parts are found.

        Returns:
            list: Geometry tags.
            Integer return if only_one was set to true.
        """
        dimtags = gmsh.model.getEntitiesInBoundingBox(x[0] - eps, y[0] - eps, z[0] - eps,
                                                      x[1] + eps, y[1] + eps, z[1] + eps,
                                                      self.dim)
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.geo_ids]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(f'Found {len(tags_filtered)} instead of only one part.')
            else:
                return tags_filtered[0]
        else:
            return tags_filtered
  
    @property
    def top_boundary(self):
        [x_min, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box([x_min, x_max], [y_max, y_max], [z_min, z_max],
                                          one_only=True)

    @property
    def bottom_boundary(self):
        [x_min, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box([x_min, x_max], [y_min, y_min], [z_min, z_max],
                                          one_only=True)

    @property
    def left_boundary(self):
        [x_min, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box([x_min, x_min], [y_min, y_max], [z_min, z_max],
                                          one_only=True)

    @property
    def right_boundary(self):
        [x_min, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box([x_max, x_max], [y_min, y_max], [z_min, z_max],
                                          one_only=True)


    def set_characteristic_length(self, char_length):
        """Set caracteristic length recursively on all boundaries and
        their boundaries.

        Args:
            char_length (float): Characteristic length for the mesh
                                 generation
        """     
        boundary = gmsh.model.getBoundary(self.dimtags, False, False, True)
        gmsh.model.mesh.setSize(boundary, char_length)

    def _make_physical(self):
        """Convert shape into physical group."""
        self.ph_id = gmsh_utils.add_physical_group(self.dim, self.geo_ids, self.name)


class MeshControl:
    """Base class for mesh restrictions."""
    
    def __init__(self, model):
        self._field = -1
        self._restricted_field = -1
        self.faces_list = []
        self.edges_list = []
        model.mesh_restrictions.append(self)

    def restrict_to_shapes(self, shapes):
        for shape in shapes:
            if shape.dim == 2:
                self.faces_list += shape.geo_ids
                self.edges_list += shape.boundaries
            if shape.dim == 1:
                self.edges_list += shape.geo_ids

    def restrict_to_faces(self, faces):
        self.faces_list += faces

    def restrict_to_edges(self, edges):
        self.edges_list += edges

    def restrict(self, shapes, faces, edges):
        self.restrict_to_shapes(shapes)
        self.restrict_to_faces(faces)
        self.restrict_to_edges(edges)

    @property
    def field(self):
        if self.faces_list == [] and self.edges_list == []:
            if self._field == -1:
                raise GmshError('Field not set.')
            return self._field
        else:
            self._restricted_field = field.add('Restrict')
            field.setNumber(self._restricted_field, 'IField', self._field)
            field.setNumbers(self._restricted_field, 'FacesList', self.faces_list)
            field.setNumbers(self._restricted_field, 'EdgesList', self.edges_list)
            return self._restricted_field


class MeshControlConstant(MeshControl):
    def __init__(self, model, char_length, shapes=[], surfaces=[], edges=[]):
        super().__init__(model)
        self._field = field.add('MathEval')
        field.setString(self._field, 'F', str(char_length))
        self.restrict(shapes, surfaces, edges)

class MeshControlLinear(MeshControl):
    def __init__(self, model, shape, min_char_length, max_char_length, dist_start=0,
                  dist_end=None, NNodesByEdge=1000, shapes=[], surfaces=[], edges=[]):
        super().__init__(model)

        if dist_end is None:
            dist_end = min_char_length + max_char_length

        if shape.dim == 1:
            edg = shape.geo_ids            
        elif shape.dim == 2:
            edg = shape.boundaries
        else:
            raise GmshError('This is only possible for shapes of dimension 1 or 2.')

        dist_field = field.add('Distance')
        field.setNumber(dist_field, 'NNodesByEdge', NNodesByEdge)
        field.setNumbers(dist_field, 'EdgesList', edg)
        self._field = field.add('Threshold')
        field.setNumber(self._field, 'IField', dist_field)
        field.setNumber(self._field, 'LcMin', min_char_length)
        field.setNumber(self._field, 'LcMax', max_char_length)
        field.setNumber(self._field, 'DistMin', dist_start)
        field.setNumber(self._field, 'DistMax', dist_end)
        self.restrict(shapes, surfaces, edges)


class MeshControlExponential(MeshControl):
    def __init__(self, model, shape, char_length, exp=1.8, fact=1, NNodesByEdge=1000,
                 shapes=[], surfaces=[], edges=[]):
        super().__init__(model)

        if shape.dim == 1:
            edg = shape.geo_ids            
        elif shape.dim == 2:
            edg = shape.boundaries
        else:
            raise GmshError('This is only possible for shapes of dimension 1 or 2.')

        dist_field = gmsh.model.mesh.field.add('Distance')
        field.setNumber(dist_field, 'NNodesByEdge', NNodesByEdge)
        field.setNumbers(dist_field, 'EdgesList', edg)
        self._field = field.add('MathEval')
        field.setString(self._field, 'F', f'F{dist_field}^{exp}*{fact}+{char_length}')
        self.restrict(shapes, surfaces, edges)
