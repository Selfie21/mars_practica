#!/usr/bin/python

import cagd.scene_2d as scene
from cagd.vec import vec2
from cagd.spline import spline, knots
from cagd.polyline import polyline
from math import pi, sin, cos
import cagd.utils as utils

#returns a list of num_samples points that are uniformly distributed on the unit circle
def unit_circle_points(num_samples):
    circle_points = []
    angle = (2 * pi) / num_samples
    for i in range(num_samples):
        circle_points.append(vec2(round(cos(i * angle), 14), round(sin(i * angle), 14)))
    return circle_points

#calculates the deviation between the given spline and a unit circle
def calculate_circle_deviation(spline):
    amount_of_edges = len(spline.control_points) - 3 # amount of control points is 3 more because 3 points are added to ensure the periodicity
    points = unit_circle_points(amount_of_edges) # points on the circle
    control_base = utils.distance(spline.control_points[0], spline.control_points[1])
    knots_base = utils.distance(points[0], points[1])
    control_area = 0.5 * control_base * amount_of_edges
    knots_area = 0.5 * knots_base * amount_of_edges
    spline_area = knots_area + ((control_area - knots_area) / 2)
    deviation = pi - spline_area
    return deviation

#interpolate 6 points with a periodic spline to create the number "8"
pts = [vec2( 0, 2.5), vec2(-1, 1), vec2( 1,-1), vec2( 0,-2.5), vec2(-1,-1), vec2(1,1)]
s = spline(3)
s = spline.interpolate_cubic_periodic(pts) 
p = s.get_polyline_from_control_points()
p.set_color("blue")
sc = scene.scene()
sc.set_resolution(900)
sc.add_element(s)
sc.add_element(p)

#generate a spline that approximates the unit circle
n = 8
circle_pts = unit_circle_points(n)
circle = spline.interpolate_cubic_periodic(circle_pts)
sc.add_element(circle)

sc.write_image()
sc.show()
