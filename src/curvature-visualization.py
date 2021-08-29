#!/usr/bin/python

#!/usr/bin/python

from cagd.spline import spline_surface, knots
from cagd.vec import vec3
import cagd.scene_2d as scene_2d
from cagd.spline import spline
from cagd.vec import vec2

pts = [vec2(0.05, 6),
       vec2(0.2, 6),
       vec2(1, 5),
       vec2(.25, 4),
       vec2(1, 3.75),
       vec2(.55, 3.5),
       vec2(.5, 3.4),
       vec2(.5, .6),
       vec2(.55, .5),
       vec2(0.8, 0.2),
       vec2(.95, 0.1),
       vec2(1, 0)]

spl = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, pts)
# you can activate these lines to view the input spline
spl.set_color("#0000ff")
sc = scene_2d.scene()
sc.set_resolution(900)
sc.add_element(spl)
sc.write_image()
# sc.show()

surface = spl.generate_rotation_surface(24)

bezier_patches = surface.to_bezier_patches()
bezier_patches.refine(1)  # refine surface into more patches for more detailed coloring
bezier_patches.visualize_curvature(bezier_patches.CURVATURE_AVERAGE, bezier_patches.COLOR_MAP_CLASSIFICATION)

path = "surfaces.off"
f = open(path, 'w')
f.write(bezier_patches.export_off())
f.close()
