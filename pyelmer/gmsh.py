"""Functions and objects that simplify the workflow with gmsh."""

from copy import deepcopy
import gmsh
import numpy as np

factory = gmsh.model.occ
field = gmsh.model.mesh.field


def cylinder(x, y, z, r, h, dim):
    """create cylinder (in 3D) or rectangle (in 2D) and return tag.

    Args:
        x (float): origin coordinate
        y (float): origin coordinate
        z (float): origin coordinate
        r (float): radius, in x-direction
        h (float): height, in y-direction
        dim (int): dimension

    Returns:
        int: tag
    """
    if dim == 2:
        return factory.addRectangle(x, y, z, r, h)
    if dim == 3:
        return factory.addCylinder(x, y, z, 0, h, 0, r)


def rotate(tag, angle=2 * np.pi, recombine=False):
    """Extrusion by rotation around y-axis. Returns tag.

    Args:
        tag (int): input surface tag
        angle (float): angle for revolution. Defaults to 2*pi.
        recombine (bool): structured mesh. Defaults to False.

    Returns:
        int: tag of 3d object
    """
    dimtags = factory.revolve(
        [(2, tag)], 0, 0, 0, 0, 1, 0, angle, heights=[0.1], recombine=recombine
    )
    for dimtag in dimtags:
        if dimtag[0] == 3:
            return dimtag[1]
    return 0  # if nothing was found; this should never happen.


def get_boundaries(dim, tag):
    """Get boundaries of gmsh entity.

    Args:
        dim (int): dimension of entity
        tag (int): tag of entity

    Returns:
        list[int]: boundary tags
    """
    dim_tags = gmsh.model.getBoundary([(dim, tag)], False, False, False)
    boundaries = []
    for dim_tag in dim_tags:
        boundaries.append(dim_tag[1])
    return boundaries


def add_physical_group(dim, tag_list, name):
    """Add physical group and set name.

    Args:
        dim (int): dimension
        tag_list (list[int]): geometry tags
        name (str): name to set

    Returns:
        int: geometry id
    """
    tag = gmsh.model.addPhysicalGroup(dim, tag_list)
    gmsh.model.setPhysicalName(dim, tag, name)
    return tag


def get_cylinder_boundary(dim, body_tag, r, h, y0, eps=1e-6):
    """Get boundary of a cylinder aligned with x-axis.
    Most useful in axi-symmetric cases.

    Args:
        dim (int): dimension
        body_tag (int): tag of entity to which cylinder belongs
        r (float): cylinder radius
        h (float): cylinder height
        y0 (float): y-coordinate of cylinder origin
        eps (float, optional): Sensitivity. Defaults to 1e-6.

    Raises:
        ValueError: if no or more than one boundary was found

    Returns:
        Int: Boundary tag
    """
    dimtags = gmsh.model.getEntitiesInBoundingBox(
        -r - eps, y0 - eps, -r - eps, r + eps, y0 + h + eps, r + eps, dim - 1
    )
    boundaries = get_boundaries(dim, body_tag)
    dimtags_filtered = [dimtag for dimtag in dimtags if dimtag[1] in boundaries]
    tags = []
    for dimtag in dimtags_filtered:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
            dimtag[0], dimtag[1]
        )
        dy = np.round(ymax - ymin, 6)
        if dy != 0:
            if np.round(xmax - r, 6) == 0:
                tags.append(dimtag[1])
    if len(tags) != 1:
        raise ValueError("Problem finding cylinder boundary :(")
    return tags[0]


def get_ring_boundary(dim, body_tag, r_out, y0, eps=1e-6):
    """Get boundary of a circle or ring around y-axis (with extent in
    x-z-plane).

    Args:
        dim (int): dimension
        body_tag (int): tag of entity to which circle / ring belongs.
        r_out (float): outer radius
        y0 (float): y-coordinate
        eps (float, optional): Sensitivity. Defaults to 1e-6.

    Raises:
        ValueError: If no or more than one boundary was found.

    Returns:
        Int: Boundary tag
    """
    dimtags = gmsh.model.getEntitiesInBoundingBox(
        -r_out - eps,
        y0 - eps,
        -r_out - eps,
        r_out + eps,
        y0 + eps,
        r_out + eps,
        dim - 1,
    )
    boundaries = get_boundaries(dim, body_tag)
    dimtags_filtered = [dimtag for dimtag in dimtags if dimtag[1] in boundaries]
    tags = []
    for dimtag in dimtags_filtered:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
            dimtag[0], dimtag[1]
        )
        dy = np.round(ymax - ymin, 6)
        if dy == 0:
            if np.round(xmax - r_out, 6) == 0:
                tags.append(dimtag[1])
    if len(tags) != 1:
        raise ValueError("Problem finding ring boundary :(")
    return tags[0]


def get_boundaries_in_box(
    x_min, y_min, z_min, x_max, y_max, z_max, dim, tag, multiple=False, eps=1e-6
):
    """Find boundaries belonging to a given entity by searching in a
    box.

    Args:
        x_min (float): min x-coordinate
        y_min (float): min y-coordinate
        z_min (float): min z-coordinate
        x_max (float): max x-coordinate
        y_max (float): max y-coordinate
        z_max (float): max z-coordinate
        dim (int): dimension of parent-entity
        tag (int): tag of parent-entity
        multiple (bool, optional): Searching for more than one boundary.
        Defaults to False.
        eps (Float, optional): Sensitivity. Defaults to 1e-6.

    Raises:
        ValueError: If more than one boundary was found (multiple=False)
        ValueError: If no boundary was found

    Returns:
        Int or list[int]: boundary tag
    """
    tags = gmsh.model.getEntitiesInBoundingBox(
        x_min - eps,
        y_min - eps,
        z_min - eps,
        x_max + eps,
        y_max + eps,
        z_max + eps,
        dim - 1,
    )
    boundaries = get_boundaries(dim, tag)
    tags_filtered = []
    for tag in tags:
        if tag[1] in boundaries:
            tags_filtered.append(tag)
    if len(tags_filtered) > 1 and not multiple:
        raise ValueError("Found more than one element!")
    elif len(tags_filtered) == 0:
        raise ValueError("Nothing found :(")
    if multiple:
        return tags_filtered  # TODO inconsistent return value
    else:
        return tags_filtered[0][1]


def set_lc(tag, lc, dim=2):
    """Set characteristic length.

    Args:
        tag (int): Tag of entity for which lc is to by set
        lc (float): Characteristic length
        dim (int, optional): Dimension. Defaults to 2.
    """
    boundary = gmsh.model.getBoundary([(dim, tag)], False, False, True)
    gmsh.model.mesh.setSize(boundary, lc)


def threshold_field(tag, lc_min, lc_max, min_dist=-1, max_dist=-1, NNodesByEdge=1000):
    """Add threshold field for mesh size control. Works in 2D only.

    Args:
        tag (int): Threshold field is applied to boundaries of the
        entity with this tag
        lc_min (float): Minimum characteristic length (applies at
        boundary)
        lc_max (float): Maximum characteristic length (far away from
        boundary)
        min_dist (int, optional): Distance to boundary where linear
        increase in element size starts. Defaults to -1.
        max_dist (int, optional): Distance to boundary where maximum
        characteristic length is reached. Defaults to -1.
        NNodesByEdge (int, optional): Resolution. Defaults to 1000.

    Returns:
        Int: Tag of threshold field.
    """
    if min_dist == -1:
        min_dist = 0
    if max_dist == -1:
        max_dist = min_dist + lc_max
    dist_field = field.add("Distance")
    field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
    field.setNumbers(dist_field, "EdgesList", get_boundaries(2, tag))
    threshold_field = field.add("Threshold")
    field.setNumber(threshold_field, "IField", dist_field)
    field.setNumber(threshold_field, "LcMin", lc_min)
    field.setNumber(threshold_field, "LcMax", lc_max)
    field.setNumber(threshold_field, "DistMin", min_dist)
    field.setNumber(threshold_field, "DistMax", max_dist)
    return threshold_field


def exp_field(tag, lc, exp=1.8, fact=1, boundaries=[], NNodesByEdge=1000):
    """Add exponential field for mesh size control. Works in 2D only.

    Args:
        tag (int): Exponential field is applied to boundaries of the
        entity with this tag.
        lc (in): Characteristic length on boundary.
        exp (float, optional): Exponent with which characteristic length
        increases. Defaults to 1.8.
        fact (float, optional): Factor in front of exponential term.
        Defaults to 1.
        boundaries (list, optional): Tags of boundaries to which
        exp-field is applied. Defaults to [].
        NNodesByEdge (int, optional): Resolution. Defaults to 1000.

    Returns:
        Int: Tag of threshold field.
    """
    if boundaries == []:
        boundaries = get_boundaries(2, tag)
    dist_field = gmsh.model.mesh.field.add("Distance")
    field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
    field.setNumbers(dist_field, "EdgesList", boundaries)
    math_field = field.add("MathEval")
    math_str = (
        "F" + str(dist_field) + "^" + str(exp) + "*" + str(fact) + " + " + str(lc)
    )
    field.setString(math_field, "F", math_str)
    return math_field


def restricted_field(base_field, faces_list=[], edges_list=[]):
    """Restrict filed to certrain faces and/or edges. Works in 2D only.

    Args:
        base_field (int): Tag of the field to be restricted
        faces_list (list, optional): Tags of faces. Defaults to [].
        edges_list (list, optional): Tags of lines. Defaults to [].

    Returns:
        Int: Tag of restricted field.
    """
    restrict_field = field.add("Restrict")
    field.setNumber(restrict_field, "IField", base_field)
    field.setNumbers(restrict_field, "FacesList", faces_list)
    field.setNumbers(restrict_field, "EdgesList", edges_list)
    return restrict_field


def restricted_const_field(surf_tag, lc, NNodesByEdge=1000):
    """Add constant field that is restricted to a certain surface.
    Works in 2D only.

    Args:
        surf_tag (int): Tag of surface.
        lc (float): Characteristic length
        NNodesByEdge (int, optional): Resolution. Defaults to 1000.

    Returns:
        Int: Tag of the resticted constant field.
    """
    math_field = field.add("MathEval")
    field.setString(math_field, "F", str(lc))
    restrict_field = field.add("Restrict")
    field.setNumber(restrict_field, "IField", math_field)
    field.setNumbers(restrict_field, "FacesList", [surf_tag])
    field.setNumbers(restrict_field, "EdgesList", get_boundaries(2, surf_tag))
    return restrict_field


def cut(obj_dimtags, tool_dimtags, remove_tool=True):
    factory.synchronize()
    out = factory.cut(obj_dimtags, tool_dimtags, removeTool=remove_tool)
    return [dimtag[1] for dimtag in out[0]]


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

    def __init__(self, name="model"):
        """Create gmsh model.

        Args:
            name (str, optional): Name of the model. Defaults to 'model'.
        """
        self._shapes = []
        self.mesh_restrictions = []
        self.min_field = -1
        self._physical = False
        gmsh.initialize(
            ["-noenv"]
        )  # see https://gitlab.onelab.info/gmsh/gmsh/-/issues/1142 for details about -noenv option
        gmsh.model.add(name)

    # TODO allow with / close -> implement __enter__, __exit__
    def close_gmsh(self):
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
        raise GeometryError(f"Shape {name} does not exist.")

    def __repr__(self):
        shapes = [s.name for s in self._shapes]
        return f"Gmsh model created with pyelmer.\nShapes: {shapes}"

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
        self.min_field = field.add("Min")
        field.setNumbers(
            self.min_field, "FieldsList", [x.field for x in self.mesh_restrictions]
        )
        field.setAsBackgroundMesh(self.min_field)

    def deactivate_characteristic_length(self):
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)

    def set_characteristic_length(self, char_length, dimensions=[0]):
        for dim in dimensions:
            gmsh.model.mesh.setSize(gmsh.model.getEntities(dim), char_length)

    def generate_mesh(self, dimension=2, order=1, size_factor=1, smoothing=1):
        self._apply_restrictions()
        gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", size_factor)
        gmsh.option.setNumber("Mesh.Smoothing", smoothing)
        gmsh.model.mesh.generate(dimension)
        gmsh.model.mesh.setOrder(order)

    def get_shapes(self, dimension, name=""):
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
            raise GmshError("This model is already physical.")
        for shape in self._shapes:
            shape._make_physical()

    def synchronize(self):
        """Synchronize gmsh geometry kernel."""
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

    def set_const_mesh_sizes(self):
        for shape in self._shapes:
            if shape.mesh_size == 0:
                print(
                    f"Warning: Mesh size = 0 for {shape.name}. Ignoring this shape..."
                )
            else:
                print(shape.name, shape.mesh_size)
                MeshControlConstant(self, shape.mesh_size, shapes=[shape])

    @property
    def symmetry_axis(self):
        sym_ax = []
        for shape in self._shapes:
            if shape.dim == 2:
                sym_ax += shape.get_boundaries_in_box([0, 0], [-1e6, 1e6])
        return sym_ax


class Shape:
    """Wrapper for any kind of shape, that shall be part of the final
    model. Shapes may be 2 or 3D objects, lines or points.
    A shape may consist of multiple gmsh entities of the same dimension.
    """

    # TODO point select of geometries / boundaries / ...?

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
        self.geo_ids = deepcopy(geo_ids)
        self.params = Parameters()
        self.ph_id = -1
        self.mesh_size = 0
        self.model._add_shape(self)

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
    def geo_id(self):
        """Gmsh geometry id.

        Raises:
            GeometryError: If there are more / less than exactly one id
            assigned to this shape.

        Returns:
            int: Gmsh tag.
        """
        if len(self.geo_ids) == 1:
            return self.geo_ids[0]
        else:
            raise GeometryError(f"This shape has {len(self.geo_ids)} geo ids.")

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
            bndry += get_boundaries(self.dim, geo_id)
        return [
            x for x in bndry if bndry.count(x) == 1
        ]  # exclude boundaries inside shape

    @property
    def bounding_box(self):
        """Get the bounding box of this shape.

        Returns:
            list[float]: [x_min, y_min, z_min, x_max, y_max, z_max]
        """
        boxes = np.array(
            [factory.getBoundingBox(self.dim, tag) for tag in self.geo_ids]
        )
        return [
            boxes[:, 0].min(),
            boxes[:, 1].min(),
            boxes[:, 2].min(),
            boxes[:, 3].max(),
            boxes[:, 4].max(),
            boxes[:, 5].max(),
        ]

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
        dimtags = gmsh.model.getEntitiesInBoundingBox(
            x[0] - eps,
            y[0] - eps,
            z[0] - eps,
            x[1] + eps,
            y[1] + eps,
            z[1] + eps,
            self.dim - 1,
        )
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.boundaries]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(
                    f"Found {len(tags_filtered)} instead of only one boundary."
                )
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
        dimtags = gmsh.model.getEntitiesInBoundingBox(
            x[0] - eps,
            y[0] - eps,
            z[0] - eps,
            x[1] + eps,
            y[1] + eps,
            z[1] + eps,
            self.dim,
        )
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.geo_ids]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(
                    f"Found {len(tags_filtered)} instead of only one part."
                )
            else:
                return tags_filtered[0]
        else:
            return tags_filtered

    @property
    def top_boundary(self):
        [x_min, _, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_max], [y_max, y_max], [z_min, z_max], one_only=True
        )

    @property
    def bottom_boundary(self):
        [x_min, y_min, z_min, x_max, _, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_max], [y_min, y_min], [z_min, z_max], one_only=True
        )

    @property
    def left_boundary(self):
        [x_min, y_min, z_min, _, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_min], [y_min, y_max], [z_min, z_max], one_only=True
        )

    @property
    def right_boundary(self):
        [_, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_max, x_max], [y_min, y_max], [z_min, z_max], one_only=True
        )

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
        self.ph_id = add_physical_group(self.dim, self.geo_ids, self.name)


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
                raise GmshError("Field not set.")
            return self._field
        else:
            self._restricted_field = field.add("Restrict")
            field.setNumber(self._restricted_field, "IField", self._field)
            field.setNumbers(self._restricted_field, "FacesList", self.faces_list)
            field.setNumbers(self._restricted_field, "EdgesList", self.edges_list)
            return self._restricted_field


class MeshControlConstant(MeshControl):
    def __init__(self, model, char_length, shapes=[], surfaces=[], edges=[]):
        super().__init__(model)
        self._field = field.add("MathEval")
        field.setString(self._field, "F", str(char_length))
        self.restrict(shapes, surfaces, edges)


class MeshControlLinear(MeshControl):
    def __init__(
        self,
        model,
        shape,
        min_char_length,
        max_char_length,
        dist_start=0,
        dist_end=None,
        NNodesByEdge=1000,
        shapes=[],
        surfaces=[],
        edges=[],
    ):
        super().__init__(model)

        if dist_end is None:
            dist_end = min_char_length + max_char_length

        if shape.dim == 1:
            edg = shape.geo_ids
        elif shape.dim == 2:
            edg = shape.boundaries
        else:
            raise GmshError("This is only possible for shapes of dimension 1 or 2.")

        dist_field = field.add("Distance")
        field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
        field.setNumbers(dist_field, "EdgesList", edg)
        self._field = field.add("Threshold")
        field.setNumber(self._field, "IField", dist_field)
        field.setNumber(self._field, "LcMin", min_char_length)
        field.setNumber(self._field, "LcMax", max_char_length)
        field.setNumber(self._field, "DistMin", dist_start)
        field.setNumber(self._field, "DistMax", dist_end)
        self.restrict(shapes, surfaces, edges)


class MeshControlExponential(MeshControl):
    def __init__(
        self,
        model,
        shape,
        char_length,
        exp=1.8,
        fact=1,
        NNodesByEdge=1000,
        shapes=[],
        surfaces=[],
        edges=[],
    ):
        super().__init__(model)

        if shape.dim == 1:
            edg = shape.geo_ids
        elif shape.dim == 2:
            edg = shape.boundaries
        else:
            raise GmshError("This is only possible for shapes of dimension 1 or 2.")

        dist_field = gmsh.model.mesh.field.add("Distance")
        field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
        field.setNumbers(dist_field, "EdgesList", edg)
        self._field = field.add("MathEval")
        field.setString(self._field, "F", f"F{dist_field}^{exp}*{fact}+{char_length}")
        self.restrict(shapes, surfaces, edges)
