"""Element shape functions, used for post processing."""

import numpy as np
from dataclasses import dataclass


@dataclass
class Node:
    """Node, contains its coordinates and Temperature."""
    x: float
    y: float
    z: float = 0
    T: float = 0

    @property
    def X(self):
        """Coordinate vector"""
        return np.array([self.x, self.y])


class Triangle:
    """Base class for triangle elements.
    """
    def __init__(self, nodes):
        self.n1 = nodes[0]
        self.n2 = nodes[1]
        self.n3 = nodes[2]

    @property
    def A(self):
        """Surface of the triangle"""
        return 0.5 * ((self.n2.x * self.n3.y - self.n3.x * self.n2.y)
                      - (self.n1.x * self.n3.y - self.n3.x * self.n1.y)
                      + (self.n1.x * self.n2.y - self.n2.x * self.n1.y))


class Triangle1st(Triangle):
    """Triangle first order.
    Node numbering:
    n3          
    |`\ 
    |  `\ 
    |    `\ 
    |      `\ 
    |        `\ 
    n1---------n2
    """
    def B_e(self, x, y):
        """Derivative of shape functions.

        Args:
            x (float): coordinate (unused for elements of this order)
            y (float): coordinate (unused for elements of this order)

        Returns:
            numpy.array: derivative matrix
        """
        del x, y
        # according to Ottosen1992 p.124 eqn 7.99
        B_e = np.array([[self.n2.y - self.n3.y, self.n3.y - self.n1.y, self.n1.y - self.n2.y],
                        [self.n3.x - self.n2.x, self.n1.x - self.n3.x, self.n2.x - self.n1.x]])
        B_e /= 2 * self.A
        return B_e

    @property
    def T(self):
        """Temperature vector"""
        T = [self.n1.T, self.n2.T, self.n3.T]
        return np.array(T)


class Triangle2nd(Triangle):
    """Triangle second order.
    Node numbering:
    n3          
    |`\ 
    |  `\ 
    n6   n5
    |      `\ 
    |        `\ 
    n1---n4----n2

    Formulas according to Zienkiewicz2005 p.116 ff.
    """
    def __init__(self, nodes):
        self.n1 = nodes[0]
        self.n2 = nodes[1]
        self.n3 = nodes[2]
        self.n4 = nodes[3]
        self.n5 = nodes[4]
        self.n6 = nodes[5]

    @property
    def a1(self):
        return self.n2.x * self.n3.y - self.n3.x * self.n2.y

    @property
    def a2(self):
        return self.n3.x * self.n1.y - self.n1.x * self.n3.y

    @property
    def a3(self):
        return self.n1.x * self.n2.y - self.n2.x * self.n1.y

    @property
    def b1(self):
        return self.n2.y - self.n3.y

    @property
    def b2(self):
        return self.n3.y - self.n1.y

    @property
    def b3(self):
        return self.n1.y - self.n2.y

    @property
    def c1(self):
        return self.n3.x - self.n2.x

    @property
    def c2(self):
        return self.n1.x - self.n3.x

    @property
    def c3(self):
        return self.n2.x - self.n1.x

    def L1(self, x, y):
        return (self.a1 + self.b1 * x + self.c1 * y)  / (2 * self.A)

    def L2(self, x, y):
        return (self.a2 + self.b2 * x + self.c2 * y)  / (2 * self.A)

    def L3(self, x, y):
        return (self.a3 + self.b3 * x + self.c3 * y)  / (2 * self.A)

    def N1(self, x, y):
        L1 = self.L1(x, y)
        return (2*L1 - 1) * L1

    def N2(self, x, y):
        L2 = self.L2(x, y)
        return (2*L2 - 1) * L2

    def N3(self, x, y):
        L3 = self.L3(x, y)
        return (2*L3 - 1) * L3

    def N4(self, x, y):
        L1 = self.L1(x, y)
        L2 = self.L2(x, y)
        return 4 * L1 * L2

    def N5(self, x, y):
        L2 = self.L2(x, y)
        L3 = self.L3(x, y)
        return 4 * L2 * L3

    def N6(self, x, y):
        L3 = self.L3(x, y)
        L1 = self.L1(x, y)
        return 4 * L3 * L1

    @property
    def L1_x(self):
        return self.b1 / (2 * self.A)

    @property
    def L1_y(self):
        return self.c1 / (2 * self.A)

    @property
    def L2_x(self):
        return self.b2 / (2 * self.A)

    @property
    def L2_y(self):
        return self.c2 / (2 * self.A)

    @property
    def L3_x(self):
        return self.b3 / (2 * self.A)

    @property
    def L3_y(self):
        return self.c3 / (2 * self.A)

    def N1_x(self, x, y):
        L1 = self.L1(x, y)
        L1_x = self.L1_x
        return 4 * L1 * L1_x - L1_x

    def N1_y(self, x, y):
        L1 = self.L1(x, y)
        L1_y = self.L1_y
        return 4 * L1 * L1_y - L1_y

    def N2_x(self, x, y):
        L2 = self.L2(x, y)
        L2_x = self.L2_x
        return 4 * L2 * L2_x - L2_x

    def N2_y(self, x, y):
        L2 = self.L2(x, y)
        L2_y = self.L2_y
        return 4 * L2 * L2_y - L2_y

    def N3_x(self, x, y):
        L3 = self.L3(x, y)
        L3_x = self.L3_x
        return 4 * L3 * L3_x - L3_x

    def N3_y(self, x, y):
        L3 = self.L3(x, y)
        L3_y = self.L3_y
        return 4 * L3 * L3_y - L3_y

    def N4_x(self, x, y):
        L1 = self.L1(x, y)
        L1_x = self.L1_x
        L2 = self.L2(x, y)
        L2_x = self.L2_x
        return 4 * (L1_x * L2 + L1 * L2_x)

    def N4_y(self, x, y):
        L1 = self.L1(x, y)
        L1_y = self.L1_y
        L2 = self.L2(x, y)
        L2_y = self.L2_y
        return 4 * (L1_y * L2 + L1 * L2_y)

    def N5_x(self, x, y):
        L2 = self.L2(x, y)
        L2_x = self.L2_x
        L3 = self.L3(x, y)
        L3_x = self.L3_x
        return 4 * (L2_x * L3 + L2 * L3_x)

    def N5_y(self, x, y):
        L2 = self.L2(x, y)
        L2_y = self.L2_y
        L3 = self.L3(x, y)
        L3_y = self.L3_y
        return 4 * (L2_y * L3 + L2 * L3_y)

    def N6_x(self, x, y):
        L3 = self.L3(x, y)
        L3_x = self.L3_x
        L1 = self.L1(x, y)
        L1_x = self.L1_x
        return 4 * (L3_x * L1 + L3 * L1_x)

    def N6_y(self, x, y):
        L3 = self.L3(x, y)
        L3_y = self.L3_y
        L1 = self.L1(x, y)
        L1_y = self.L1_y
        return 4 * (L3_y * L1 + L3 * L1_y)

    def B_e(self, x, y):
        """Derivative of shape functions.

        Args:
            x (float): coordinate
            y (float): coordinate

        Returns:
            numpy array: Derivative matrix at (x, y)
        """
        B_e = [
            [self.N1_x(x, y), self.N2_x(x, y), self.N3_x(x, y), self.N4_x(x, y), self.N5_x(x, y), self.N6_x(x, y)],
            [self.N1_y(x, y), self.N2_y(x, y), self.N3_y(x, y), self.N4_y(x, y), self.N5_y(x, y), self.N6_y(x, y)]
        ]
        return np.array(B_e)

    @property
    def T(self):
        """Temperature vector"""
        T = [self.n1.T, self.n2.T, self.n3.T, self.n4.T, self.n5.T, self.n6.T]
        return np.array(T)

    @T.setter
    def T(self, T):
        """Temperature vector"""
        self.n1.T = T[0]
        self.n2.T = T[1]
        self.n3.T = T[2]
        self.n4.T = T[3]
        self.n5.T = T[4]
        self.n6.T = T[5]


class Line1st:
    """Line element of first order:
    n1----n2
    """
    def __init__(self, nodes):
        self.n1 = nodes[0]
        self.n2 = nodes[1]
        n = np.array([self.n2.y - self.n1.y, -(self.n2.x - self.n1.x)])
        self.normal = n / np.linalg.norm(n)
    
    def invert_normal(self):
        self.normal *= -1


class Line2nd(Line1st):
    """Line element of second order:
    n1----n3----n2
    """
    def __init__(self, nodes):
        super().__init__(nodes[:-1])
        self.n3 = nodes[2]
