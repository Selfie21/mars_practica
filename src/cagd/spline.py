#! /usr/bin/python
import copy
from math import *

import cagd.utils as utils
from cagd.bezier import bezier_patches
from cagd.polyline import polyline
from cagd.vec import vec2, vec3


class spline:
    # Interpolation modes
    INTERPOLATION_EQUIDISTANT = 0
    INTERPOLATION_CHORDAL = 1
    INTERPOLATION_CENTRIPETAL = 2
    INTERPOLATION_FOLEY = 3

    def __init__(self, degree):
        assert (degree >= 1)
        self.degree = degree
        self.knots = None
        self.control_points = []
        self.color = "black"

    # checks if the number of knots, controlpoints and degree define a valid spline
    def validate(self):
        knots = self.knots.validate()
        points = len(self.knots) == len(self.control_points) + self.degree + 1
        return knots and points

    def evaluate(self, t):
        a, b = self.support()
        assert (a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            # the spline is only defined on the interval [a, b)
            # it is useful to define self(b) as lim t->b self(t)
            t = t - 0.000001
        return self.de_boor(t, 1)[0]

    # returns the interval [a, b) on which the spline is supported
    def support(self):
        return self.knots[self.degree], self.knots[len(self.knots) - self.degree - 1]

    def __call__(self, t):
        return self.evaluate(t)

    def tangent(self, t):
        a, b = self.support()
        assert (a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            t = t - 0.000001
        return self.de_boor(t, 2)[1] - self.de_boor(t, 2)[0]

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    # calculates the de_boor scheme at a given value u
    # stops when the column is only "stop" elements long
    # returns that column as a list
    def de_boor(self, u, stop):
        n = self.degree
        t = self.knots
        i = self.knots.knot_index(u)
        j = i - n

        offset = n + 1
        len_points = n
        d_length = int(((n + 1) * (n + 2)) / 2)
        assert i >= n and i + n < len(self.knots)
        d = [0] * d_length
        for tmp in range(0, offset):
            d[tmp] = self.control_points[j + tmp]

        # using one array to save storage, using offset to keep track of current position in array
        # tmpj is the index for the d in the kth iteration (0-len_points)
        for k in range(1, n + 1):
            for tmpj in range(len_points):
                j = i - n + tmpj
                alpha = (u - t[j + k]) / (t[j + n + 1] - t[j + k])
                d[tmpj + offset] = (1 - alpha) * d[tmpj + offset - len_points - 1] + (
                        alpha * d[tmpj + offset - len_points])
            offset += len_points
            len_points -= 1

            if len_points < stop:
                return d[offset - (len_points + 1):offset]

    # adjusts the control points such that it represents the same function,
    # but with an added knot
    def insert_knot(self, t):
        index = self.knots.knot_index(t)
        prev_controlpoints = copy.deepcopy(self.control_points)
        for i in range(index - 2, index):
            alpha = (t - self.knots[i]) / (self.knots[i + 3] - self.knots[i])
            self.control_points[i] = (1 - alpha) * prev_controlpoints[i - 1] + alpha * prev_controlpoints[i]

        alpha = (t - self.knots[index]) / (self.knots[index + 3] - self.knots[index])
        d = (1 - alpha) * prev_controlpoints[index - 1] + alpha * prev_controlpoints[index]
        self.control_points.insert(index, d)
        self.knots.insert(t)

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
            # print("comparing {0} to {1} and {2}".format(p, min_vec, max_vec))
            if p.x < min_vec.x:
                min_vec.x = p.x
            if p.y < min_vec.y:
                min_vec.y = p.y
            if p.x > max_vec.x:
                max_vec.x = p.x
            if p.y > max_vec.y:
                max_vec.y = p.y
        return min_vec, max_vec

    def draw(self, scene, num_samples):
        i = self.degree - 1
        while i < len(self.knots) - self.degree - 2:
            i += 1
            k0 = self.knots[i]
            k1 = self.knots[i + 1]
            if k0 == k1:
                continue
            p0 = self(k0)
            for j in range(1, num_samples + 1):
                t = k0 + j / num_samples * (k1 - k0)
                p1 = self(t)
                scene.draw_line(p0, p1, self.color)
                p0 = p1

    def get_polyline_from_control_points(self):
        pl = polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl

    # generates a spline that interpolates the given points using the given mode
    # returns that spline object
    def interpolate_cubic(mode, points):
        new_spline = spline(3)
        new_spline.initialize_knots(points)
        new_spline.generate_knots(mode, points)
        diag1, diag2, diag3, resx, resy = new_spline.generate_sole(points)
        control_points_x = utils.solve_tridiagonal_equation(diag1, diag2, diag3, resx)
        control_points_y = utils.solve_tridiagonal_equation(diag1, diag2, diag3, resy)
        for pt_x, pt_y in zip(control_points_x, control_points_y):
            new_spline.control_points.append(vec2(pt_x, pt_y))
        return new_spline

    # Calculates and returns αi
    def alpha(self, i):
        return (self.knots[i + 2] - self.knots[i]) / (self.knots[i + 3] - self.knots[i])

    # Calculates and returns βi
    def beta(self, i):
        return (self.knots[i + 2] - self.knots[i + 1]) / (self.knots[i + 3] - self.knots[i + 1])

    # Calculates γi
    def gamma(self, i):
        return (self.knots[i + 2] - self.knots[i + 1]) / (self.knots[i + 4] - self.knots[i + 1])

    # generates the system of linear equations for creating a spline
    def generate_sole(self, points):
        dim = len(points) + 2
        diag1 = [0] * dim
        diag2 = [0] * dim
        diag3 = [0] * dim

        diag1[1] = -1
        for i in range(2, dim - 2):
            diag1[i] = (1 - self.beta(i)) * (1 - self.alpha(i))
        diag1[-2] = -1 + self.gamma(dim - 3)
        diag1[-1] = 0

        diag2[0] = 1
        diag2[1] = 1 + self.alpha(2)
        for i in range(2, dim - 2):
            diag2[i] = (1 - self.beta(i)) * self.alpha(i) + self.beta(i) * (1 - self.gamma(i))
        diag2[-2] = -self.gamma(dim - 3) + 2
        diag2[-1] = 1

        diag3[0] = 0
        diag3[1] = -self.alpha(2)
        for i in range(2, dim - 2):
            diag3[i] = self.beta(i) * self.gamma(i)
        diag3[-2] = -1

        resx, resy = self.generate_residuum(dim - 2, points)
        resx.insert(1, 0)
        resx.insert(-1, 0)
        resy.insert(1, 0)
        resy.insert(-1, 0)
        return diag1, diag2, diag3, resx, resy

    def generate_residuum(self, dim, points):
        resx = [0] * dim
        resy = [0] * dim
        for i in range(dim):
            resx[i] = points[i].x

        for i in range(dim):
            resy[i] = points[i].y

        return resx, resy

    def generate_knots(self, mode, points):
        m = len(self.knots)

        if mode == spline.INTERPOLATION_EQUIDISTANT:
            for i in range(1, m):
                self.knots[i] = i
        elif mode == spline.INTERPOLATION_CHORDAL:
            for i in range(1, m):
                prev_point = points[i - 1]
                current_point = points[i]
                self.knots[i] = utils.distance(prev_point, current_point) + self.knots[i - 1]
        elif mode == spline.INTERPOLATION_CENTRIPETAL:
            for i in range(1, m):
                prev_point = points[i - 1]
                current_point = points[i]
                self.knots[i] = sqrt(utils.distance(prev_point, current_point)) + self.knots[i - 1]
        elif mode == spline.INTERPOLATION_FOLEY:
            for i in range(0, m - 1):
                prev_point = points[i - 1]
                current_point = points[i]
                next_point = points[i + 1]

                theta = min(pi - utils.angle(prev_point, current_point), pi / 2)
                theta_next = min(pi - utils.angle(current_point, next_point), pi / 2)

                prev_d = 0 if i == 1 else utils.distance(points[i - 1], prev_point)
                current_d = utils.distance(prev_point, current_point)
                next_d = 0 if i == m - 1 else utils.distance(current_point, next_point)
                self.knots[i + 1] = (current_d * (1
                                                  + ((3 * theta * prev_d) / ((2 * prev_d) + current_d))
                                                  + ((3 * theta_next * next_d) / ((2 * next_d) + current_d)))) + \
                                    self.knots[i]
        self.quadruple_edge_knots()

    # intializes knots with first and last point
    def initialize_knots(self, points):
        self.knots = knots(len(points))
        self.knots[0] = 0

    # quadruple edge knots
    def quadruple_edge_knots(self):
        for i in range(3):
            self.knots.insert(self.knots[0])
            self.knots.insert(self.knots[-1])

    # generates a spline that interpolates the given points and fulfills the definition
    # of a periodic spline
    # returns that spline object
    def interpolate_cubic_periodic(points):
        new_spline = spline(3)
        points.append(points[0])
        dim = len(points) - 1

        new_spline.knots = knots(dim + 7)
        new_spline.knots[0] = 0
        for i in range(1, dim + 7):
            new_spline.knots[i] = i

        diag1 = [(1 / 6)] * dim
        diag2 = [(4 / 6)] * dim
        diag3 = [(1 / 6)] * dim

        resx, resy = new_spline.generate_residuum(dim, points)
        control_points_x = utils.solve_almost_tridiagonal_equation(diag1, diag2, diag3, resx)
        control_points_y = utils.solve_almost_tridiagonal_equation(diag1, diag2, diag3, resy)
        for pt_x, pt_y in zip(control_points_x, control_points_y):
            new_spline.control_points.append(vec2(pt_x, pt_y))

        new_spline.convert_to_periodic_controlpoints(control_points_x, control_points_y)
        return new_spline

    def convert_to_periodic_controlpoints(self, control_points_x, control_points_y):
        self.control_points.insert(0, vec2(control_points_x[len(control_points_x) - 1],
                                           control_points_y[len(control_points_y) - 1]))
        self.control_points.append(vec2(control_points_x[0], control_points_y[0]))
        self.control_points.append(vec2(control_points_x[1], control_points_y[1]))
        self.control_points.append(vec2(control_points_x[2], control_points_y[2]))
        self.control_points.append(vec2(control_points_x[3], control_points_y[3]))

    # for splines of degree 3, generate a parallel spline with distance dist
    # the returned spline is off from the exact parallel by at most eps
    def generate_parallel(self, dist, eps):
        assert (self.degree == 3)
        para_spline = self.interpolate_parallel_spline(dist)

        while not self.distance_threshold_between_splines(dist, eps, para_spline):
            para_spline = self.interpolate_parallel_spline(dist)

        return para_spline

    # returns of the distance inbetween two knot points between this spline and the parallel spline is smaller than eps
    # if it isn't adds a knot in between these two knots
    def distance_threshold_between_splines(self, distance, eps, parallel_spline):
        for i in range(3, len(self.knots) - 3):
            inbetween = self.knots[i] + (self.knots[i + 1] - self.knots[i]) / 2
            normal_pt = self.evaluate(inbetween)
            para_pt = parallel_spline.evaluate(inbetween)
            current_distance = utils.distance(normal_pt, para_pt)
            if abs((current_distance - abs(distance))) > eps:
                self.insert_knot(inbetween)
                return False
        return True

    # creates a parallel spline
    def interpolate_parallel_spline(self, dist):
        interpolation_points = [self.evaluate(knot) for knot in self.knots[3:-3]]
        index = 0

        for knot in self.knots[3:-3]:
            tangent = self.tangent(knot)
            normal = vec2(tangent.y, -tangent.x)
            normal *= 1 / utils.euklidian_norm(normal)
            interpolation_points[index] += (dist * normal)
            index += 1

        para_spline = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, interpolation_points)
        para_spline.knots = copy.deepcopy(self.knots)
        return para_spline

    # generates a rotational surface by rotating the spline around the z axis
    # the spline is assumed to be on the xz-plane
    # num_samples refers to the number of interpolation points in the rotational direction
    # returns a spline surface object in three dimensions
    def generate_rotation_surface(self, num_samples):

        surface = spline_surface((2, 3))

        for ctrl_point in self.control_points:

            rotation_points = []
            for i in range(num_samples):
                rotation_point = self.generate_rotation_point(ctrl_point, i, num_samples)
                rotation_points.append(rotation_point)

            control_points_periodic3d = []
            circle_spline = spline.interpolate_cubic_periodic(rotation_points)
            control_points_periodic2d = circle_spline.control_points
            z = ctrl_point.y
            for control_point2d in control_points_periodic2d:
                net_point = vec3(control_point2d.x, control_point2d.y, z)
                control_points_periodic3d.append(net_point)

            surface.knots = (self.knots, circle_spline.knots)
            surface.control_points.append(control_points_periodic3d)
        return surface

    def generate_rotation_point(self, ctrl_point, index, num_samples):
        return vec3(ctrl_point.x * cos((2 * pi * index) / num_samples),
                    ctrl_point.x * sin((2 * pi * index) / num_samples),
                    ctrl_point.y)


class spline_surface:
    # the two directions of the parameter space
    DIR_U = 0
    DIR_V = 1

    # creates a spline of degrees m,n
    # degree is a tuple (m,n)
    def __init__(self, degree):
        du, dv = degree
        assert (du >= 1 and dv >= 1)
        self.degree = degree
        self.knots = (None, None)  # tuple of both knot vectors
        self.control_points = []  # 2dim array of control points

    # checks if the number of knots, controlpoints and degree define a valid spline
    def validate(self):
        if len(self.control_points) == 0:
            return False
        k1, k2 = self.knots
        d1, d2 = self.degree
        knots = k1.validate() and k2.validate()
        p1 = len(self.control_points)
        p2 = len(self.control_points[0])
        points1 = len(k1) == p1 + d1 + 1
        points2 = len(k2) == p2 + d2 + 1
        return knots and points1 and points2

    def evaluate(self, u, v):
        s1, s2 = self.support()
        a, b = s1
        c, d = s2
        assert (a <= u <= b and c <= v <= v)
        if u == b:
            u = u - 0.000001
        if v == d:
            v = v - 0.000001
        t = (u, v)
        return self.de_boor(t, (1, 1))[0][0]

    # return nested tuple ((a,b), (c,d))
    # the spline is supported in (u,v) \in [a,b)x[c,d]
    def support(self):
        k1, k2 = self.knots
        d1, d2 = self.degree
        s1 = (k1[d1], k1[len(k1) - d1 - 1])
        s2 = (k2[d2], k2[len(k2) - d2 - 1])
        return (s1, s2)

    def __call__(self, u, v):
        return self.evaluate(u, v)

    # calculates the de boor scheme at t = (u,v)
    # until there are only stop = (s1, s2) elements left
    def de_boor(self, t, stop):
        d1, d2 = self.degree
        k1, k2 = self.knots
        s1, s2 = stop
        u, v = t
        m1 = len(self.control_points)
        m2 = len(self.control_points[0])

        new_rows = [None for i in range(m1)]
        for row in range(m1):
            spl = spline(d2)
            spl.knots = k2
            spl.control_points = self.control_points[row]
            new_rows[row] = spl.de_boor(v, s2)

        new_pts = [None for i in range(s2)]
        for col in range(s2):
            spl = spline(d1)
            spl.knots = k1
            ctrl_pts = [new_rows[i][col] for i in range(m1)]
            spl.control_points = ctrl_pts
            new_pts[col] = spl.de_boor(u, s1)

        return new_pts

    def insert_knot(self, direction, t):
        if direction == self.DIR_U:
            self.__insert_knot_u(t)
        elif direction == self.DIR_V:
            self.__insert_knot_v(t)
        else:
            assert (False)

    def __insert_knot_v(self, t):
        du, dv = self.degree
        ku, kv = self.knots
        nu = len(self.control_points)
        nv = len(self.control_points[0])
        for i in range(nu):
            row = self.control_points[i]
            spl = spline(du)
            spl.control_points = copy.copy(row)
            spl.knots = copy.deepcopy(kv)
            spl.insert_knot(t)
            self.control_points[i] = spl.control_points
        kv.insert(t)

    def __insert_knot_u(self, t):
        du, dv = self.degree
        ku, kv = self.knots
        nu = len(self.control_points)
        nv = len(self.control_points[0])
        new_control_points = [[None for i in range(nv)] for j in range(nu + 1)]
        for i in range(nv):
            col = [self.control_points[j][i] for j in range(nu)]
            spl = spline(dv)
            spl.control_points = col
            spl.knots = copy.deepcopy(ku)
            spl.insert_knot(t)
            for j in range(nu + 1):
                new_control_points[j][i] = spl.control_points[j]
        self.control_points = new_control_points
        ku.insert(t)

    def to_bezier_patches(self):
        patches = bezier_patches()
        #m, n = 3
        m=3
        n=3
        for i in range(len(self.knots[0])):
            p=0
            if(self.knots[0][i] == self.knots[0][i+1]):
                p = p + 1
            else:
                for j in range(m - p):
                    self.insert_knot(self.DIR_U, self.knots[0][i])
                i = i + 1
        for i in range(len(self.knots[1])):
            p=0
            if(self.knots[1][i] == self.knots[1][i+1]):
                p = p + 1
            else:
                for j in range(n - p):
                    self.insert_knot(self.DIR_V, self.knots[1][i])
                i = i +p
        return patches


class knots:
    # creates a knots array with n elements
    def __init__(self, n):
        self.knots = [None for i in range(n)]

    def validate(self):
        prev = None
        for k in self.knots:
            if k is None:
                return False
            if prev is None:
                prev = k
            else:
                if k < prev:
                    return False
        return True

    def __len__(self):
        return len(self.knots)

    def __getitem__(self, i):
        return self.knots[i]

    def __setitem__(self, i, v):
        self.knots[i] = v

    def __delitem__(self, i):
        del self.knots[i]

    def __iter__(self):
        return iter(self.knots)

    def insert(self, t):
        i = 0
        while self[i] < t:
            i += 1
        self.knots.insert(i, t)

    def knot_index(self, v):
        if v < self.knots[0] or v >= self.knots[-1]:
            return -1

        for index in range(self.__len__() - 1):
            if v < self.knots[index + 1]:
                return index
