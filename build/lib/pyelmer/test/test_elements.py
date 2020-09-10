# created by Arved Enders-Seidlitz on 07.09.2020
#
# tests for elements.py

import copy
from pyelmer.elements import Node, Triangle2nd, Line2nd

# simple input data
n1 = Node(0, 0, 0, 0)
n2 = Node(2, 0, 0, 4)
n3 = Node(0, 2, 0, 1)
n4 = Node(1, 0, 0, 1)
n5 = Node(1, 1, 0, 1)
n6 = Node(0, 1, 0, 1)
tr = Triangle2nd([n1, n2, n3, n4, n5, n6])

def values_at_nodes(N, nodes):
    values = []
    for n in nodes:
        values.append(N(n.x, n.y))
    return values

def test_shape_functions():
    assert tr.N1(n1.x, n1.y) == 1
    assert tr.N2(n2.x, n2.y) == 1
    assert tr.N3(n3.x, n3.y) == 1
    assert tr.N4(n4.x, n4.y) == 1
    assert tr.N5(n5.x, n5.y) == 1
    assert tr.N6(n6.x, n6.y) == 1

    assert values_at_nodes(tr.N1, [n2, n3, n4, n5, n6]) == [0, 0, 0, 0, 0]
    assert values_at_nodes(tr.N2, [n1, n3, n4, n5, n6]) == [0, 0, 0, 0, 0]
    assert values_at_nodes(tr.N3, [n1, n2, n4, n5, n6]) == [0, 0, 0, 0, 0]
    assert values_at_nodes(tr.N4, [n1, n2, n3, n5, n6]) == [0, 0, 0, 0, 0]
    assert values_at_nodes(tr.N5, [n1, n2, n3, n4, n6]) == [0, 0, 0, 0, 0]
    assert values_at_nodes(tr.N6, [n1, n2, n3, n4, n5]) == [0, 0, 0, 0, 0]

def test_derivatives():
    assert tr.N1_x(n1.x, n1.y) == tr.N1_y(n1.x, n1.y)
    assert tr.N2_x(n2.x, n2.x) == tr.N3_y(n3.x, n3.y)
    assert tr.N2_y(n2.x, n2.x) == tr.N3_x(n3.x, n3.y)
    assert tr.N4_x(n4.x, n4.y) == 0
    assert tr.N6_y(n6.x, n6.y) == 0
    assert tr.N4_y(n4.x, n4.y) == tr.N6_x(n6.x, n6.y)
    assert tr.N5_x(n5.x, n5.y) == tr.N5_y(n5.x, n5.y)

def test_line():
    n1 = Node(0, 0, 0)
    n2 = Node(1, 0, 0)
    n3 = Node(2, 0, 0)
    l = Line2nd([n1, n2, n3])
    normal = copy.deepcopy(l.normal)
    l.invert_normal()
    assert list(l.normal) == list(-1* normal)

if __name__ == "__main__":
    test_shape_functions()
    test_derivatives()
    test_line()
