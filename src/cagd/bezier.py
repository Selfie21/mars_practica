#!/usr/bin/python

import copy
from math import sqrt

import cagd.utils as utils
from cagd.polyline import polyline


class bezier_curve:
    def __init__(self, degree):
        assert (degree >= 0)
        self.degree = degree
        self.control_points = [None for i in range(degree + 1)]
        self.color = "black"

    def set_control_point(self, index, val):
        assert (index >= 0 and index <= self.degree)
        self.control_points[index] = val

    def get_control_point(self, index):
        assert (index >= 0 and index <= self.degree)
        return self.control_points[index]

    # evaluates the curve at t
    def evaluate(self, t):
        return self.__de_casteljeau(t, 1)[0]

    # evaluates tangent at t
    def tangent(self, t):
        last_two_ctrl_pts = self.__de_casteljeau(t, 2)
        a = last_two_ctrl_pts[0]
        b = last_two_ctrl_pts[1]
        return b - a

    # calculates the normal at t
    def normal(self, t):
        pass

    # syntactic sugar so bezier curve can be evaluated as curve(t)
    # instead of curve.evaluate(t)
    def __call__(self, t):
        return self.evaluate(t)

    # calculates the de-casteljeau scheme until the column only has stop elements
    def __de_casteljeau(self, t, stop):
        assert (stop >= 1)
        column = self.control_points
        while len(column) > stop:
            new_column = [None for i in range(len(column) - 1)]
            for i in range(len(new_column)):
                new_column[i] = (1 - t) * column[i] + t * column[i + 1]
            column = new_column
        return column

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    # calculates the bezier representation of the derivative
    def get_derivative(self):
        pass

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
            if p.x < min_vec.x:
                min_vec.x = p.x
            if p.y < min_vec.y:
                min_vec.y = p.y
            if p.x > max_vec.x:
                max_vec.x = p.x
            if p.y > max_vec.y:
                max_vec.y = p.y
        return (min_vec, max_vec)

    def draw(self, scene, num_samples):
        p0 = self(0)
        for i in range(1, num_samples + 1):
            t = i / num_samples
            p1 = self(t)
            scene.draw_line(p0, p1, self.color)
            p0 = p1

    def get_polyline_from_control_points(self):
        pl = polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl


class bezier_surface:
    # creates a bezier surface of degrees n,m
    # the degree parameter is a tuple (n,m)
    def __init__(self, degree):
        d1, d2 = degree
        assert (d1 >= 0 and d2 >= 0)
        self.degree = degree
        self.control_points = [[None for i in range(d2 + 1)] for j in range(d1 + 1)]
        white = (1, 1, 1)
        self.color = (white, white, white, white)
        self.curvature = (None, None, None, None)

    def set_control_point(self, index1, index2, val):
        assert (index1 >= 0 and index1 <= self.degree[0])
        assert (index2 >= 0 and index2 <= self.degree[1])
        self.control_points[index1][index2] = val

    def get_control_point(self, index1, index2):
        assert (index1 >= 0 and index1 <= self.degree[0])
        assert (index2 >= 0 and index2 <= self.degree[1])
        return self.control_points[index1][index2]

    def evaluate(self, t1, t2):
        return self.__de_casteljeau(t1, t2, (1, 1))[0][0]

    # sets the colors at the corners
    # c00 is the color at u=v=0, c01 is the color at u=0 v=1, etc
    # a color is a tuple (r,g,b) with values between 0 an 1
    def set_colors(self, c00, c01, c10, c11):
        self.color = (c00, c01, c10, c11)

    # sets the curvature at the corners
    # c00 is the curvature at u=v=0, c01 is the curvature at u=0 v=1, etc
    def set_curvature(self, c00, c01, c10, c11):
        self.curvature = (c00, c01, c10, c11)

    def get_curvature(self):
        return self.curvature

    def __call__(self, t):
        t1, t2 = t
        return self.evaluate(t1, t2)

    def __de_casteljeau(self, t1, t2, stop):
        s1, s2 = stop
        d1, d2 = self.degree
        assert (s1 >= 1 and s2 >= 1)
        d1 += 1  # number of control points in each direction
        d2 += 1

        # apply the casteljeau scheme in one direction,
        # ie, reduce dimension from (d1, d2) to (s1, d2)
        column = self.control_points
        while d1 > s1:
            d1 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t1) * column[i][j] + t1 * column[i + 1][j]
            column = new_column

        # apply the casteljeau scheme in the other direction,
        # ie, reduce dimension from (s1, d2) to (s1, s2)
        while d2 > s2:
            d2 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t2) * column[i][j] + t2 * column[i][j + 1]
            column = new_column

        return column

    def normal(self, t1, t2):
        pass

    def get_derivative(self, direction):
        pass

    def subdivide(self, t1, t2):
        b0, b1 = self.__subdivide_u(t1)
        b00, b01 = b0.__subdivide_v(t2)
        b10, b11 = b1.__subdivide_v(t2)
        return [b00, b01, b10, b11]

    def __subdivide_u(self, t):
        du, dv = self.degree
        left = bezier_surface((du, dv))
        right = bezier_surface((du, dv))
        for k in range(du + 1):
            pts = self.__de_casteljeau(t, 0, (du - k + 1, dv + 1))
            left.control_points[k] = pts[0]
            right.control_points[k] = pts[-1]
        return (left, right)

    def __subdivide_v(self, t):
        du, dv = self.degree
        left = bezier_surface((du, dv))
        right = bezier_surface((du, dv))
        for k in range(dv + 1):
            pts = self.__de_casteljeau(0, t, (du + 1, dv - k + 1))
            for i in range(dv + 1):
                left.control_points[i][k] = pts[i][0]
                right.control_points[i][k] = pts[i][-1]
        return (left, right)

    def calculate_colors(surface, color_map, min_value, max_value):
        c = 4 * [(0, 0, 0)]
        for i in range(4):
            x = surface.get_curvature()[i]
            if color_map == bezier_patches.COLOR_MAP_CUT:
                if x < 0:
                    x = 0
                elif x > 1:
                    x = 1
            elif color_map == bezier_patches.COLOR_MAP_LINEAR:
                x = (x - min_value) / (max_value - min_value)
            elif color_map == bezier_patches.COLOR_MAP_CLASSIFICATION:
                if x < 0:
                    x = 0
                elif x == 0:
                    x = 0.5
                else:
                    x = 1
            else:
                raise ValueError('Color map not supported')

            if 0 <= x <= 0.25:
                c[i] = (0, 4 * x, 1)
            elif 0.25 < x <= 0.5:
                c[i] = (0, 1, 2 - (4 * x))
            elif 0.5 < x <= 0.75:
                c[i] = ((4 * x)-2, 1, 0)
            elif 0.75 < x <= 1:
                c[i] = (1, 4 - (4 * x), 0)
            else:
                raise ValueError('Curvature value calculation went wrong')
        return c[0], c[1], c[2], c[3]


class bezier_patches:
    CURVATURE_GAUSSIAN = 0
    CURVATURE_AVERAGE = 1
    CURVATURE_PRINCIPAL_MAX = 2  # Maximale Hauptkruemmung
    CURVATURE_PRINCIPAL_MIN = 3  # Minimale Hauptkruemmung
    COLOR_MAP_LINEAR = 4
    COLOR_MAP_CUT = 5
    COLOR_MAP_CLASSIFICATION = 6

    def __init__(self):
        self.patches = []

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, p):
        return self.patches[p]

    def __setitem__(self, i, p):
        self.patches[i] = p

    def __delitem__(self, p):
        del self.patches[p]

    def __iter__(self):
        return iter(self.patches)

    def append(self, p):
        self.patches.append(p)

    # refines patches by subdividing each patch into four new patches
    # there are 4^num times more patches after calling this function
    def refine(self, num):
        for i in range(num):
            new_patches = bezier_patches()
            for p in self:
                new = p.subdivide(0.5, 0.5)
                for n in new:
                    new_patches.append(n)
            self.patches = new_patches

    def visualize_curvature(self, curvature_mode, color_map):
        all_curvature = []
        for surface in self.patches:
            current_surface_curvature = [bezier_patches.calculate_curvature(surface, curvature_mode, x, y) for x in
                                         range(2) for y in range(2)]
            surface.set_curvature(*current_surface_curvature)
            all_curvature.append(current_surface_curvature)

        min_value = min(map(min, all_curvature))
        max_value = max(map(max, all_curvature))
        for surface in self.patches:
            c00, c01, c10, c11 = surface.calculate_colors(color_map, min_value, max_value)
            surface.set_colors(c00, c01, c10, c11)

    @staticmethod
    def calculate_curvature(patch, curvature_mode, x, y):
        b_u, b_u_u, b_v, b_v_v, b_u_v = bezier_patches.get_all_derivatives(patch, x, y)
        n = b_u.cross_product(b_v) * (1 / utils.euclidean_norm(b_u.cross_product(b_v)))
        g = [b_u.dot(b_u), b_u.dot(b_v), b_v.dot(b_u), b_v.dot(b_v)]
        h = [n.dot(b_u_u), n.dot(b_u_v), n.dot(b_u_v), n.dot(b_v_v)]

        gaussian = utils.determinant_m2(h) / utils.determinant_m2(g)
        average = ((h[0] * g[3]) - (2 * h[1] * g[1]) + (h[3] * g[0])) / (2 * ((g[0] * g[3]) - (g[1] * g[1])))

        if curvature_mode == bezier_patches.CURVATURE_GAUSSIAN:
            return gaussian
        elif curvature_mode == bezier_patches.CURVATURE_AVERAGE:
            return average
        elif curvature_mode == bezier_patches.CURVATURE_PRINCIPAL_MAX:
            return average + sqrt((average * average) - gaussian)
        elif curvature_mode == bezier_patches.CURVATURE_PRINCIPAL_MIN:
            return average - sqrt((average * average) - gaussian)
        else:
            raise ValueError('Curvature Mode not available!')

    @staticmethod
    def get_all_derivatives(patch, x, y):
        b_u = bezier_patches.partial_derivatives(patch, 'u', 1)[-x][-y]
        b_u_u = bezier_patches.partial_derivatives(patch, 'u', 2)[-x][-y]
        b_v = bezier_patches.partial_derivatives(patch, 'v', 1)[-x][-y]
        b_v_v = bezier_patches.partial_derivatives(patch, 'v', 2)[-x][-y]
        b_u_v = bezier_patches.partial_derivatives(patch, 'both', 2)[-x][-y]
        return b_u, b_u_u, b_v, b_v_v, b_u_v

    # creates the k_th derivative in the corresponding direction for a specific patch
    @staticmethod
    def partial_derivatives(patch, direction, k):
        m = patch.degree[0]
        n = patch.degree[1]
        u_factor = 0
        v_factor = 0
        scale = 0
        control_points = patch.control_points

        if direction == 'u':
            u_factor = k
            scale = m
        elif direction == 'v':
            v_factor = k
            scale = n
        else:
            u_factor = int(k / 2)
            v_factor = int(k / 2)

        deriv = [[None for col in range((n + 1 - v_factor))] for row in range(m + 1 - u_factor)]
        for i in range(0, m + 1 - u_factor):
            for j in range(0, n + 1 - v_factor):
                # single case
                if u_factor + v_factor <= 1:
                    deriv[i][j] = scale * (control_points[i + u_factor][j + v_factor] - control_points[i][j])
                # multi case
                elif direction == 'u':
                    deriv[i][j] = scale * (m - 1) * (
                            control_points[i + 2][j] - (2 * control_points[i + 1][j]) + control_points[i][j])
                elif direction == 'v':
                    deriv[i][j] = scale * (n - 1) * (
                            control_points[i][j + 2] - (2 * control_points[i][j + 1]) + control_points[i][j])
                else:
                    deriv[i][j] = n * m * (
                            control_points[i + 1][j + 1] - control_points[i][j + 1] - control_points[i + 1][j] +
                            control_points[i][j])
        return deriv

    def export_off(self):
        def export_point(p):
            return str(p.x) + " " + str(p.y) + " " + str(p.z)

        def export_colors(c):
            s = ""
            for x in c:
                s += str(x)
                s += " "
            s += "1"  # opacity
            return s

        s = "CBEZ333\n"
        for patch in self:
            # coordinates
            for row in patch.control_points:
                for p in row:
                    s += export_point(p)
                    s += "\n"

            # colors
            s += export_colors(patch.color[0])
            s += "\n"
            s += export_colors(patch.color[2])
            s += "\n"
            s += export_colors(patch.color[1])
            s += "\n"
            s += export_colors(patch.color[3])
            s += "\n"
            s += "\n"

        return s
