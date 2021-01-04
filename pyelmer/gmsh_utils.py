"""Utility functions for gmsh"""

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


def rotate(tag):
    """Extrusion by rotation around y-axis. Returns tag.

    Args:
        tag (int): input surface tag

    Returns:
        int: tag of 3d object
    """
    dimtags = factory.revolve([(2, tag)], 0, 0, 0, 0, 1, 0, 2*np.pi)
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
    dimtags = gmsh.model.getEntitiesInBoundingBox(-r - eps, y0 - eps, -r - eps, r + eps,
                                                  y0 + h + eps, r + eps, dim-1)
    boundaries = get_boundaries(dim, body_tag)
    dimtags_filtered = [dimtag for dimtag in dimtags if dimtag[1] in boundaries]
    tags = []
    for dimtag in dimtags_filtered:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dimtag[0], dimtag[1])
        dy = np.round(ymax - ymin, 6)
        if dy != 0:
            if np.round(xmax - r, 6) == 0:
                tags.append(dimtag[1])
    if len(tags) != 1:
        raise ValueError('Problem finding cylinder boundary :(')
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
    dimtags = gmsh.model.getEntitiesInBoundingBox(-r_out -eps, y0 - eps, -r_out - eps,
                                                  r_out + eps, y0 + eps, r_out + eps, dim-1)
    boundaries = get_boundaries(dim, body_tag)
    dimtags_filtered = [dimtag for dimtag in dimtags if dimtag[1] in boundaries]
    tags = []
    for dimtag in dimtags_filtered:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dimtag[0], dimtag[1])
        dy = np.round(ymax - ymin, 6)
        if dy == 0:
            if np.round(xmax - r_out, 6) == 0:
                tags.append(dimtag[1])
    if len(tags) != 1:
        raise ValueError('Problem finding ring boundary :(')
    return tags[0]


def get_boundaries_in_box(x_min, y_min, z_min, x_max, y_max, z_max, dim, tag, multiple=False,
                          eps=1e-6):
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
    tags = gmsh.model.getEntitiesInBoundingBox(x_min - eps, y_min - eps, z_min - eps, x_max + eps,
                                               y_max + eps, z_max + eps, dim-1)
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
        return tags_filtered
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
    dist_field = field.add('Distance')
    field.setNumber(dist_field, 'NNodesByEdge', NNodesByEdge)
    field.setNumbers(dist_field, 'EdgesList', get_boundaries(2, tag))
    threshold_field = field.add('Threshold')
    field.setNumber(threshold_field, 'IField', dist_field)
    field.setNumber(threshold_field, 'LcMin', lc_min)
    field.setNumber(threshold_field, 'LcMax', lc_max)
    field.setNumber(threshold_field, 'DistMin', min_dist)
    field.setNumber(threshold_field, 'DistMax', max_dist)
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
    dist_field = gmsh.model.mesh.field.add('Distance')
    field.setNumber(dist_field, 'NNodesByEdge', NNodesByEdge)
    field.setNumbers(dist_field, 'EdgesList', boundaries)
    math_field = field.add('MathEval')
    math_str = 'F' + str(dist_field) + '^' + str(exp) + '*' + str(fact) + ' + ' + str(lc)
    field.setString(math_field, 'F', math_str)
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
    restrict_field = field.add('Restrict')
    field.setNumber(restrict_field, 'IField', base_field)
    field.setNumbers(restrict_field, 'FacesList', faces_list)
    field.setNumbers(restrict_field, 'EdgesList', edges_list)
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
    math_field = field.add('MathEval')
    field.setString(math_field, 'F', str(lc))
    restrict_field = field.add('Restrict')
    field.setNumber(restrict_field, 'IField', math_field)
    field.setNumbers(restrict_field, 'FacesList', [surf_tag])
    field.setNumbers(restrict_field, 'EdgesList', get_boundaries(2, surf_tag))
    return restrict_field


def cut(obj_dimtags, tool_dimtags, remove_tool=True):
    factory.synchronize()
    out = factory.cut(obj_dimtags, tool_dimtags, removeTool=remove_tool)
    tags = []
    for dimtag in out[0]:
        tags.append(dimtag[1])
    return tags
