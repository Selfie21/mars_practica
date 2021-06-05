#!/usr/bin/python

import numpy as np

import cagd.scene_2d as scene_2d
# create an example spline to demonstrate how to create a spline
# you can use this to test your implementation of the de-boor algorithm
#    and the knot_index function
from cagd.polyline import polyline
from cagd.spline import spline, knots
from cagd.utils import solve_almost_tridiagonal_equation
from cagd.utils import solve_tridiagonal_equation
from cagd.vec import vec2

example_spline = spline(3)
example_spline.control_points = [vec2(0, 0), vec2(1, 1), vec2(2, 2), vec2(3, 3), vec2(4, 4), vec2(5, 5), vec2(6, 6)]
example_spline.knots = knots(11)
example_spline.knots.knots = [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1]
example_spline.de_boor(0.4, 1)
p = example_spline.get_polyline_from_control_points()
p.set_color("red")

# test for solve tridiagonal equation operation
A = np.array([[11, 2, 0, 0, 0], [3, 9, 5, 0, 0], [0, 6, 25, 8, 0], [0, 0, 3, 26, 7], [0, 0, 0, 4, 13]], dtype=float)
a = [0, 3, 6, 3, 4]
b = [11, 9, 25, 26, 13]
c = [2, 5, 8, 7, 0]
d = [5, 7, 8, 4, 3]
print("Test results:")
print(np.around(solve_tridiagonal_equation(a, b, c, d), decimals=8))
print(np.linalg.solve(A, d))

# test for solve almost tridiagonal equation operation
A = np.array([[11, 2, 0, 0, 3], [3, 9, 5, 0, 0], [0, 6, 25, 8, 0], [0, 0, 3, 26, 7], [5, 0, 0, 4, 13]], dtype=float)
a = [3, 3, 6, 3, 4]
b = [11, 9, 25, 26, 13]
c = [2, 5, 8, 7, 5]
d = [5, 7, 8, 4, 3]
print("\nTest results 2:")
print(np.around(solve_almost_tridiagonal_equation(a, b, c, d), decimals=8))
print(np.linalg.solve(A, d))

# interpolate six points with the four different interpolation options to
#    draw a small letter "e"
# uncomment these lines once you implemented the spline interpolation
pts = [vec2(0, .4), vec2(.8, .8), vec2(.5, 1.2), vec2(-.03, .4), vec2(.4, 0), vec2(1, .2)]
s1 = spline.interpolate_cubic(spline.INTERPOLATION_EQUIDISTANT, pts)
s2 = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, pts)
s3 = spline.interpolate_cubic(spline.INTERPOLATION_CENTRIPETAL, pts)
s4 = spline.interpolate_cubic(spline.INTERPOLATION_FOLEY, pts)
s1.set_color("#000066")
s2.set_color("#0000aa")
s3.set_color("#6666ff")
s4.set_color("#aaaaff")
p = polyline()
p.points = pts
p.set_color("red")

# generate a scene and add elements to it
sc = scene_2d.scene()
sc.set_resolution(900)
sc.add_element(p)
sc.add_element(s1)
sc.add_element(s2)
sc.add_element(s3)
sc.add_element(s4)
sc.write_image()  # compose all elements in the scene
sc.show()  # tries to show the image with a default viewer
sc.write_to_file("test.png")  # saves the image to a file
