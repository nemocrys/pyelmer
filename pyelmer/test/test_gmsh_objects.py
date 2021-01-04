from dataclasses import dataclass
import pytest
import gmsh
from pyelmer.gmsh_objects import *

@dataclass
class Rects:
    r1: Shape
    r2: Shape
    r3: Shape

@pytest.fixture(scope='module')
def rectangles():
    with Model() as model:
        r1 = factory.addRectangle(0, 0, 0, 1, 1)
        r2 = factory.addRectangle(1, 0, 0, 1, 1)
        factory.fragment([(2, r1)], [(2, r2)])
        factory.synchronize()
        rect1 = Shape(model, 2, 'r1', [r1])
        rect2 = Shape(model, 2, 'r2', [r2])
        rect3 = Shape(model, 2, 'r3', [r1, r2])
        yield Rects(rect1, rect2, rect3)


def test_geo_ids(rectangles):
    assert rectangles.r1.geo_ids == [1]
    assert rectangles.r2.geo_ids == [2]
    assert rectangles.r3.geo_ids == [1, 2]


def test_get_boundaries(rectangles):
    assert sorted(rectangles.r1.get_boundaries()) == [1, 2, 3, 4]
    assert sorted(rectangles.r2.get_boundaries()) == [2, 5, 6, 7]
    assert sorted(rectangles.r3.get_boundaries()) == [1, 3, 4, 5, 6, 7]


if __name__ == "__main__":
    test_get_boundaries(rectangles)