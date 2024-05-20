import pymol
from pymol import cmd
from wills_pymol_crap.pymol_util import getres, com

CHAINS = r"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"

def cgo_arrow(xyz1, xyz2, radius=0.5, col=(1, 1, 1), gap=0.0, hlength=-1, hradius=-1, name=''):
   from chempy import cpv

   xyz1 = [xyz1.x, xyz1.y, xyz1.z]
   xyz2 = [xyz2.x, xyz2.y, xyz2.z]

   radius, gap = float(radius), float(gap)
   hlength, hradius = float(hlength), float(hradius)

   color1 = list(col)
   color2 = list(col)
   # xyz1 = get_coord(atom1)
   # xyz2 = get_coord(atom2)
   normal = cpv.normalize(cpv.sub(xyz1, xyz2))

   if hlength < 0:
      hlength = radius * 3.0
   if hradius < 0:
      hradius = hlength * 0.6

   if gap:
      diff = cpv.scale(normal, gap)
      xyz1 = cpv.sub(xyz1, diff)
      xyz2 = cpv.add(xyz2, diff)

   xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

   obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
         [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
         [1.0, 0.0]

   if not name:
      name = cmd.get_unused_name('arrow')

   return obj
   # cmd.load_cgo(obj, name)

def connectors(sele, name='connectors', frac=(0.3, 0.3), radius=1):
   cmd.delete(name)
   # chains = cmd.get_chains(sele)
   objs = cmd.get_object_list(sele)
   cgo = list()

   # coms = [com(f'({sele}) and chain {c}') for c in chains]
   coms = [com(f'({sele}) and {object}') for object in objs]
   coms = sorted(coms, key=lambda com: com.z)

   # ic(coms)
   # return

   for acom, bcom in zip(coms[:-1], coms[1:]):
      a = acom * (1 - frac[0]) + bcom * frac[0]
      b = acom * frac[1] + bcom * (1 - frac[1])
      # cgo += cgo_cyl(a, b, radius, col=(0.5,0.5,0.5))
      cgo += cgo_arrow(a, b, radius, col=(0.5, 0.5, 0.5), hlength=10, hradius=5)
   cmd.load_cgo(cgo, name)

def expandhelix(sele, raddelta=0, lenscale=1):
   for c in cmd.get_chains(sele):
      sel = f'({sele}) and (chain {c})'
      ic(sel)
      cen = com(sel)
      delta = Vec(cen.x, cen.y, 0).normalized() * raddelta
      delta.z = cen.z * (lenscale - 1)
      ic(c, delta)
      trans(sel, delta)

def alignmono(mono, template):
   # resch = getres('hbnet')
   # res = [r[1] for r in resch]
   # print(res)
   for c in cmd.get_chains(template):
      if c in '!#$@': continue
      newname = f'{mono}_{c}'
      cmd.create(newname, mono)
      # cmd.alter(newname, f'chain="{c}"')
      cmd.super(newname, f'({template}) and chain {c}')
   # sel = f'({mono}_*) and (resi {str(res)[1:-1].replace( ", ","+")})'
   # cmd.select('allhbnet',f'{sel} and not name n+ca+c+o')

def alignmono_ab(mono, template):
   # resch = getres('hbnet')
   # res = [r[1] for r in resch]
   # print(res)
   ichain = 0
   for c in cmd.get_chains(template):
      newname = f'{mono}_{c}'
      cmd.create(newname, mono)
      # cmd.alter(f'{newname} and chain A',f'chain="{CHAINS[ichain]}"')
      # ichain += 1
      # cmd.alter(f'{newname} and chain B',f'chain="{CHAINS[ichain]}"')
      # ichain += 1
      cmd.super(newname, f'({template}) and chain {c}')
   # sel = f'({mono}_*) and (resi {str(res)[1:-1].replace( ", ","+")})'
   # cmd.select('allhbnet',f'{sel} and not name n+ca+c+o')