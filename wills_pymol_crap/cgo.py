import math
import pymol
from pymol.cgo import *
import random
from wills_pymol_crap.xyzMath import *
import numpy as np

def cgo_circle(
      cen,
      r=1.0,
      w=2.0,
      xform=Xform(),
      col=(1, 1, 1),
):

   x = float(cen.x)
   y = float(cen.y)
   z = float(cen.z)
   r = abs(float(r))
   cr = abs(float(col[0]))
   cg = abs(float(col[1]))
   cb = abs(float(col[2]))
   w = float(w)

   obj = [BEGIN, LINES, COLOR, cr, cg, cb]
   for i in np.arange(360) * np.pi / 180:
      pt1 = Vec(r * math.cos(i), r * math.sin(i), 0)
      pt2 = Vec(r * math.cos(i + np.pi / 180), r * math.sin(i + np.pi / 180), 0)
      pt1 = xform * pt1 + cen
      pt2 = xform * pt2 + cen
      # print(pt1, pt2)
      obj.append(VERTEX)
      obj.append(pt1.x)
      obj.append(pt1.y)
      obj.append(pt1.z)
      obj.append(VERTEX)
      obj.append(pt2.x)
      obj.append(pt2.y)
      obj.append(pt2.z)
   obj.append(END)

   cName = cmd.get_unused_name("circle_")
   cmd.load_cgo(obj, cName)
   cmd.set("cgo_line_width", w, cName)
   return obj
