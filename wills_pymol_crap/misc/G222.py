from random import random
from pymol.cgo import *
from wills_pymol_crap.misc.gyroid_coords import points as gxyz, triangles as gtri, normals as gnorm

def gyroid(s=62.35, fx=-0.25, fy=0.25, fz=-0.25, r=80, c=Vec(0, 0, 0), t=0.1, a=1.0):
   cmd.set('suspend_updates','on')
   cmd.delete('gyroid')
   v = cmd.get_view()
   r2 = r * r
   obj = [BEGIN, TRIANGLES, ALPHA, a]
   for dx in (-1, 0, 1):
      for dy in (-1, 0, 1):
         for dz in (-1, 0, 1):
            for i, j, k in gtri:
               pi0, pi1, pi2 = s * gxyz[i][0] + s * (dx + fx), s * gxyz[i][1] + s * (
                  dy + fy), s * gxyz[i][2] + s * (dz + fz)
               pj0, pj1, pj2 = s * gxyz[j][0] + s * (dx + fx), s * gxyz[j][1] + s * (
                  dy + fy), s * gxyz[j][2] + s * (dz + fz)
               pk0, pk1, pk2 = s * gxyz[k][0] + s * (dx + fx), s * gxyz[k][1] + s * (
                  dy + fy), s * gxyz[k][2] + s * (dz + fz)
               ni0, ni1, ni2 = gnorm[i][0], gnorm[i][1], gnorm[i][2]
               nj0, nj1, nj2 = gnorm[j][0], gnorm[j][1], gnorm[j][2]
               nk0, nk1, nk2 = gnorm[k][0], gnorm[k][1], gnorm[k][2]
               pcen = Vec((pi0 + pj0 + pk0) / 3.0, (pi1 + pj1 + pk1) / 3.0,
                          (pi2 + pj2 + pk2) / 3.0)
               if c.distance_squared(pcen) > r2: continue
               obj.extend(
                  (COLOR, 1, 0, 1, NORMAL, ni0, ni1, ni2, VERTEX, pi0 + t * ni0, pi1 + t * ni1,
                   pi2 + t * ni2, NORMAL, nj0, nj1, nj2, VERTEX, pj0 + t * nj0, pj1 + t * nj1,
                   pj2 + t * nj2, NORMAL, nk0, nk1, nk2, VERTEX, pk0 + t * nk0, pk1 + t * nk1,
                   pk2 + t * nk2, COLOR, 1, 1, 0, NORMAL, -ni0, -ni1, -ni2, VERTEX, pi0 - t * ni0,
                   pi1 - t * ni1, pi2 - t * ni2, NORMAL, -nk0, -nk1, -nk2, VERTEX, pk0 - t * nk0,
                   pk1 - t * nk1, pk2 - t * nk2, NORMAL, -nj0, -nj1, -nj2, VERTEX, pj0 - t * nj0,
                   pj1 - t * nj1, pj2 - t * nj2))
   obj.append(END)
   cmd.load_cgo(obj, 'gyroid')
   cmd.set_view(v)
   cmd.set('suspend_updates','off')   

def gyroid_wireframe(s=62.35, fx=-0.25, fy=0.25, fz=-0.25, r=80, c=Vec(0, 0, 0), a=1.0):
   cmd.set('suspend_updates','on')
   v = cmd.get_view()
   r2 = r * r
   obj = [BEGIN, LINES, COLOR, 1, 1, 1, ALPHA, a]
   count = 0
   for dx in (-1.01, 0, 1.01):
      for dy in (-1.01, 0, 1.01):
         for dz in (-1.01, 0, 1.01):
            for i, j, k in gtri:
               pi0, pi1, pi2 = s * gxyz[i][0] + s * (dx + fx), s * gxyz[i][1] + s * (
                  dy + fy), s * gxyz[i][2] + s * (dz + fz)
               pj0, pj1, pj2 = s * gxyz[j][0] + s * (dx + fx), s * gxyz[j][1] + s * (
                  dy + fy), s * gxyz[j][2] + s * (dz + fz)
               pk0, pk1, pk2 = s * gxyz[k][0] + s * (dx + fx), s * gxyz[k][1] + s * (
                  dy + fy), s * gxyz[k][2] + s * (dz + fz)
               pcen = Vec((pi0 + pj0 + pk0) / 3.0, (pi1 + pj1 + pk1) / 3.0,
                          (pi2 + pj2 + pk2) / 3.0)
               if c.distance_squared(pcen) > r2: continue
               obj.extend((COLOR, (dx + 2.0) / 3.0, (dy + 2.0) / 3.0, (dz + 2.0) / 3.0, VERTEX,
                           pi0, pi1, pi2, VERTEX, pj0, pj1, pj2, VERTEX, pj0, pj1, pj2, VERTEX,
                           pk0, pk1, pk2, VERTEX, pk0, pk1, pk2, VERTEX, pi0, pi1, pi2))
            count += 1
   obj.append(END)
   cmd.load_cgo(obj, 'gyroid')
   cmd.set_view(v)
   cmd.set('suspend_updates','off')   
