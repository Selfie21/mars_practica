#!/usr/bin/python

from cagd.polyline import polyline
from cagd.spline import spline
from cagd.vec import vec2
import cagd.scene_2d as scene_2d

pts = [vec2(0,.4), vec2(.8,.8), vec2(.5,1.2), vec2(-.03,.4), vec2(.4,0), vec2(1,.2)]
s1 = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, pts)
para = s1.generate_parallel(0.1, 0.005)
s1.set_color("#0000ff")

sc = scene_2d.scene()
sc.set_resolution(900)
sc.add_element(para)
sc.add_element(s1)
sc.write_image()
sc.show()
