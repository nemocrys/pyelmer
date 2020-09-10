# created by Arved Enders-Seidlitz on 14.07.2020
#
# utilities for gmsh

import gmsh
import numpy as np
factory = gmsh.model.occ

def get_boundary(dim, tag, coordinate, min_max, boundaries=[], sym=False, filter_y=False):
    """find the boundary with minimum / maximum location in selected coordinate

    Args:
        dim (int): dimension of element
        tag (int): tag of element
        coordinate (int): 0 - x coordinate, 1 - y coordinate, 2 - z coordinate
        min_max (str): 'min' or 'max' - minimum or maximum in coordinate.
        boundaries (list): use these boundaries instead of those given by tag

    Returns:
        int: tag of boundary element or 0 in case of an error
    """
    if min_max != 'min' and min_max != 'max':
        raise ValueError("Select between 'min' or 'max'")
    if boundaries == []:
        boundaries = gmsh.model.getBoundary([dim, tag], oriented=False)
    coord_min = []
    coord_max = []
    dy = []
    for boundary in boundaries:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(boundary[0], boundary[1])
        x_min = np.round([xmin, ymin, zmin], 6)
        x_max = np.round([xmax, ymax, zmax], 6)
        coord_min.append(x_min[coordinate])
        coord_max.append(x_max[coordinate])
        dy.append(x_max[1] - x_min[1])
    if min_max == 'min':
        if sym:
            if filter_y:
                for i in range(len(coord_min)):
                    if dy[i] == 0:
                        coord_min[i] = 1e6
            i = [i for i, x in enumerate(coord_min) if x == min(coord_min)]
        else:
            i = [i for i, x in enumerate(coord_max) if x == min(coord_max)]
        if len(i) != 1:
            return 0
        return boundaries[i[0]][1]
    else:
        if dim == 3 and sym:
            if filter_y:
                for i in range(len(coord_max)):
                    if dy[i] == 0:
                        coord_max[i] = -1e6
            i = [i for i, x in enumerate(coord_max) if x == max(coord_max)]
        else:
            i = [i for i, x in enumerate(coord_min) if x == max(coord_min)]
        if len(i) != 1:
            return 0
        return boundaries[i[0]][1]

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
    dim_tags = gmsh.model.getBoundary([(dim, tag)], False, False, False)
    boundaries = []
    for dim_tag in dim_tags:
        boundaries.append(dim_tag[1])
    return boundaries

def add_physical_group(dim, tag_list, name):
    tag = gmsh.model.addPhysicalGroup(dim, tag_list)
    gmsh.model.setPhysicalName(dim, tag, name)
    return tag

def get_cylinder_boundary(dim, body_tag, r, h, y0, eps=1e-6):
    dimtags = gmsh.model.getEntitiesInBoundingBox(-r - eps, y0 - eps, -r - eps,
        r + eps, y0 + h + eps, r + eps, dim-1)
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


        
def get_element_in_box(x_min, y_min, z_min, x_max, y_max, z_max, dim, tag, multiple=False, eps = 1e-6):
    tags = gmsh.model.getEntitiesInBoundingBox(x_min - eps, y_min - eps, z_min - eps,
        x_max + eps, y_max + eps, z_max + eps, dim-1)
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
    boundary = gmsh.model.getBoundary([(dim, tag)], False, False, True)
    gmsh.model.mesh.setSize(boundary, lc)

def threshold_field(tag, lc_min, lc_max, min_dist=-1, max_dist=-1, NNodesByEdge=1000):
    # function works in 2D only
    field = gmsh.model.mesh.field
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
    field = gmsh.model.mesh.field
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
    field = gmsh.model.mesh.field
    restrict_field = field.add('Restrict')
    field.setNumber(restrict_field, 'IField', base_field)
    field.setNumbers(restrict_field, 'FacesList', faces_list)
    field.setNumbers(restrict_field, 'EdgesList', edges_list)
    return restrict_field

def restricted_const_field(surf_tag, lc, NNodesByEdge=1000):
    field = gmsh.model.mesh.field
    math_field = field.add('MathEval')
    field.setString(math_field, 'F', str(lc))
    restrict_field = field.add('Restrict')
    field.setNumber(restrict_field, 'IField', math_field)
    field.setNumbers(restrict_field, 'FacesList', [surf_tag])
    field.setNumbers(restrict_field, 'EdgesList', get_boundaries(2, surf_tag))
    return restrict_field

