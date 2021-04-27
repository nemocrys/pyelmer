from dataclasses import dataclass
import pytest
import gmsh
from pyelmer.gmsh import *


@dataclass
class Rects:
    r1: Shape
    r2: Shape
    r3: Shape


@pytest.fixture(scope="module")
def rectangles():
    model = Model()
    r1 = factory.addRectangle(0, 0, 0, 1, 1)
    r2 = factory.addRectangle(1, 0, 0, 1, 1)
    factory.fragment([(2, r1)], [(2, r2)])
    factory.synchronize()
    rect1 = Shape(model, 2, "r1", [r1])
    rect2 = Shape(model, 2, "r2", [r2])
    rect3 = Shape(model, 2, "r3", [r1, r2])
    # gmsh.fltk.run()
    yield Rects(rect1, rect2, rect3)


def test_geo_ids(rectangles):
    assert rectangles.r1.geo_ids == [1]
    assert rectangles.r2.geo_ids == [2]
    assert rectangles.r3.geo_ids == [1, 2]


def test_get_boundaries(rectangles):
    assert sorted(rectangles.r1.boundaries) == [1, 2, 3, 4]
    assert sorted(rectangles.r2.boundaries) == [2, 5, 6, 7]
    assert sorted(rectangles.r3.boundaries) == [1, 3, 4, 5, 6, 7]


def test_get_boundaries_in_box(rectangles):
    assert sorted(rectangles.r1.get_boundaries_in_box([-0.5, 1.5], [-0.5, 1.5])) == [
        1,
        2,
        3,
        4,
    ]
    assert rectangles.r1.get_boundaries_in_box([0.5, 1.5], [-0.5, 1.5]) == [2]
    assert (
        rectangles.r1.get_boundaries_in_box([0.5, 1.5], [-0.5, 1.5], one_only=True) == 2
    )


def test_top_bottom_left_right_boundary(rectangles):
    assert rectangles.r1.top_boundary == 3
    assert rectangles.r1.bottom_boundary == 1
    assert rectangles.r1.right_boundary == 2
    assert rectangles.r1.left_boundary == 4

    assert rectangles.r2.top_boundary == 7
    assert rectangles.r2.bottom_boundary == 5
    assert rectangles.r2.right_boundary == 6
    assert rectangles.r2.left_boundary == 2
