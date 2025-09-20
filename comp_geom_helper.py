import math
import numpy as np

class TetrahedronCalculations:
    """
    Provides geometric calculations for a tetrahedron defined by four 3D points.
    """

    def __init__(self, points):
        """
        Initialize with an array-like of four 3D points.
        """
        self.points = np.asarray(points)

    def tetra_edge_lengths(self):
        """
        Returns the lengths of all edges, their sum, and average.
        """
        a, b, c, d = self.points
        lengths = [
            np.linalg.norm(a - b),
            np.linalg.norm(a - c),
            np.linalg.norm(a - d),
            np.linalg.norm(b - c),
            np.linalg.norm(b - d),
            np.linalg.norm(c - d)
        ]
        sum_len = sum(lengths)
        avg_len = sum_len / 6
        return (*lengths, sum_len, avg_len)

    def tetra_volume(self):
        """
        Returns the volume of the tetrahedron using the scalar triple product.
        """
        a, b, c, d = self.points
        return abs(np.dot(a - d, np.cross(b - d, c - d))) / 6

    def triangle_area(self, p1, p2, p3):
        """
        Returns the area of a triangle given three 3D points.
        """
        return 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))

    def tetra_surface_area(self):
        """
        Returns the areas of the four faces and the total surface area.
        """
        a, b, c, d = self.points
        t1 = self.triangle_area(a, b, c)
        t2 = self.triangle_area(a, b, d)
        t3 = self.triangle_area(a, c, d)
        t4 = self.triangle_area(b, c, d)
        tetra_area = t1 + t2 + t3 + t4
        return t1, t2, t3, t4, tetra_area

    def solid_angle(self, v, p1, p2, p3):
        """
        Returns the solid angle at vertex v subtended by triangle (p1, p2, p3).
        """
        a = p1 - v
        b = p2 - v
        c = p3 - v
        numerator = abs(np.dot(a, np.cross(b, c)))
        denom = (
            np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c)
            + np.dot(a, b) * np.linalg.norm(c)
            + np.dot(a, c) * np.linalg.norm(b)
            + np.dot(b, c) * np.linalg.norm(a)
        )
        return 2 * np.arctan2(numerator, denom)

    def tetra_solid_angles(self):
        """
        Returns the solid angles at each vertex of the tetrahedron.
        """
        a, b, c, d = self.points
        return (
            self.solid_angle(a, b, c, d),
            self.solid_angle(b, a, c, d),
            self.solid_angle(c, a, b, d),
            self.solid_angle(d, a, b, c)
        )

class TriangleCalculations:
    """
    Provides geometric calculations for a triangle defined by three 3D points.
    """

    def __init__(self, points):
        """
        Initialize with an array-like of three 3D points.
        """
        self.points = np.asarray(points)

    def distance(self, p1, p2):
        """
        Returns the Euclidean distance between two 3D points.
        """
        return np.linalg.norm(p1 - p2)

    def triangle_properties(self):
        """
        Returns the side lengths and area of the triangle.
        """
        a = self.distance(self.points[0], self.points[1])
        b = self.distance(self.points[1], self.points[2])
        c = self.distance(self.points[2], self.points[0])
        s = (a + b + c) / 2
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        return a, b, c, area

    @staticmethod
    def triangle_area(p1, p2, p3):
        """
        Returns the area of a triangle given three 3D points.
        """
        return 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))

def consecutive_among_three(nums):
    """
    Returns the length of the longest sequence of consecutive numbers among three numbers.
    """
    def max_consecutive_adjacent(lst):
        max_count = 0
        current_count = 1
        for i in range(1, len(lst)):
            if lst[i] - lst[i - 1] == 1:
                current_count += 1
                max_count = max(max_count, current_count)
            else:
                current_count = 1
        return max_count if max_count >= 2 else 0
    return max(max_consecutive_adjacent(nums), max_consecutive_adjacent(nums[::-1]))

def consecutive_among_four(nums):
    """
    Returns the length of the longest sequence of consecutive numbers among four numbers.
    """
    def max_consecutive_adjacent(lst):
        max_count = 0
        current_count = 1
        for i in range(1, len(lst)):
            if lst[i] - lst[i - 1] == 1:
                current_count += 1
                max_count = max(max_count, current_count)
            else:
                current_count = 1
        return max_count if max_count >= 2 else 0
    return max(max_consecutive_adjacent(nums), max_consecutive_adjacent(nums[::-1]))

def get_flat_list(nested_list):
    """
    Flattens a nested list into a single list.
    """
    return [item for sublist in nested_list for item in sublist]