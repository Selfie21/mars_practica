#!/usr/bin/python

import cagd.scene_2d as scene_2d
from cagd.spline import spline
from cagd.vec import vec2

pts = [vec2(0, .4), vec2(.8, .8), vec2(.5, 1.2), vec2(-.03, .4), vec2(.4, 0), vec2(1, .2)]
s1 = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, pts)
para = s1.generate_parallel(0.1, 0.005)
s1.set_color("#0000ff")

sc = scene_2d.scene()
sc.set_resolution(900)
sc.add_element(para)
sc.add_element(s1)

for i in [-1, 1]:
    para = s1.generate_parallel(i * 0.025, 0.005)
    para.set_color("#999999")
    sc.add_element(para)

sc.write_image()
sc.show()

# KIT
ptsK = [vec2(0, 1), vec2(0, 0), vec2(.25, .5), vec2(.5, 1), vec2(.25, .5), vec2(.5, 0)]
ptsI = [vec2(0, 1), vec2(0.35, 1), vec2(.1725, 1), vec2(.1725, .5), vec2(.1725, 0), vec2(0, 0), vec2(.35, 0)]
ptsI = [x + vec2(.75, 0) for x in ptsI]
ptsT = [vec2(.5, 0), vec2(.5, .5), vec2(.5, 1), vec2(1, 1), vec2(0, 1)]
ptsT = [x + vec2(1.25, 0) for x in ptsT]

sK = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, ptsK)
paraK = spline.generate_parallel(sK, 0.1, 0.005)
sI = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, ptsI)
sT = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, ptsT)
paraK.set_color("#009682")

sc = scene_2d.scene()
sc.set_resolution(900)
sc.add_element(sK)
sc.add_element(sI)
sc.add_element(sT)
sc.add_element(paraK)
sc.write_image()
sc.show()
