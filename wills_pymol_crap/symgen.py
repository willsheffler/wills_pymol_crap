import itertools, functools, re, os, inspect
from wills_pymol_crap.xyzMath import (Vec, Mat, Xform, RAD, projperp, SYMTET, SYMOCT, SYMICS, isvec, randnorm, Ux, Uy,
                                      Uz)
from wills_pymol_crap.pymol_util import ROSETTA_CHAINS, xform, rot, trans
from wills_pymol_crap.pymol_util import cgo_cyl, cgo_sphere, pymol, pdb_format_atom
from pymol import cmd, util
from wills_pymol_crap.symgen_classes import (
   SymElem,
   SymElemPosition,
   SymElemGroupManager,
   SymTrieNode,
   SymTrieSanityCheckVisitor,
   BuildCGO,
   ComponentCenterVisitor,
   VecDict,
   RosettaSymDef,
)


def get_chain_frames(sel='visible'):
   xforms = list()
   for x in cmd.get_object_list(sel):
      for ch in cmd.get_chains(x):
         f = getframe(f'{sel} and chain {ch}')
         x = np.eye(4)
         x[0,:3] = f.R.rowx()
         x[1,:3] = f.R.rowy()
         x[2,:3] = f.R.rowz()                        
         x[:3,3] = f.t
         xforms.append(x)
   xforms = np.stack(xforms)
   ic(xforms.shape)
   np.save('T4.npy', xforms)

def make_helix_simple(n=20, nstates=1):
   cmd.delete('MAKESYM')
   cmd.set('suspend_updates', 'on')
   v = cmd.get_view()

   for istate in range(max(1, nstates)):
      state = istate + 1

      frame0 = getframe('chain A and vis', state=state)
      if state in (8, 17, 19, 20):
         frame1 = getframe('chain D and vis', state=state)
      else:
         frame1 = getframe('chain B and vis', state=state)
      fdelta = frame1 * ~frame0
      frames = [Xform()]
      for i in range(n):
         frames.append(fdelta * frames[-1])

      srcstate = state if nstates > 0 else 0
      tgtstate = -1 if nstates > 0 else 0
      makesym(frames, f'vis and chain A', srcstate=srcstate, tgtstate=tgtstate)

   cmd.set_view(v)
   cmd.set('suspend_updates', 'off')

def make_p632(*args, depth=6, **kw):
   make_layer(3, 2, 'ABC', 'AD', *args, depth=depth, **kw)

def make_p333(*args, depth=4, **kw):
   make_layer(3, 3, 'ABC', 'ADE', *args, depth=depth, **kw)

def make_p442(*args, depth=5, **kw):
   make_layer(4, 2, 'ABCD', 'AE', *args, depth=depth, **kw)

def make_p444(*args, depth=4, **kw):
   make_layer(4, 4, 'ABCD', 'AEFG', *args, depth=depth, **kw)

def make_layer(
   nfold1,
   nfold2,
   cnfold1,
   cnfold2,
   sele=None,
   asym='vis and chain A',
   depth=6,
   maxrad=9e9,
   nstates=0,
   **kw,
):
   cmd.set('suspend_updates', 'on')
   v = cmd.get_view()

   cmd.delete('obj MAKESYM')
   cmd.delete('obj TMP_*')
   if sele is None:
      sele = cmd.get_object_list()[0]

   for istate in range(max(1, nstates)):
      state = istate + 1
      cen_nf1 = com(f'({sele}) and chain {"+".join(cnfold1)}', state=state)
      cen_nf2 = com(f'({sele}) and chain {"+".join(cnfold2)}', state=state)
      # generators1 = [
      # # rotation_around_degrees(Vec(0, 0, 1), 0, cen_nf1),
      # rotation_around_degrees(Vec(0, 0, 1), 120, cen_nf1),
      # rotation_around_degrees(Vec(0, 0, 1), 240, cen_nf1),
      # ]
      angs1 = [np.pi * 2 / nfold1 * i for i in range(1, nfold1)]
      angs2 = [np.pi * 2 / nfold2 * i for i in range(1, nfold2)]
      generators1 = [rotation_around(Vec(0, 0, 1), a, cen_nf1) for a in angs1]
      generators2 = [rotation_around(Vec(0, 0, 1), a, cen_nf2) for a in angs2]
      # generators2 = [
      # # rotation_around_degrees(Vec(0, 0, 1), 0, cen_nf2),
      # rotation_around_degrees(Vec(0, 0, 1), 180, cen_nf2),
      # ]
      g = generators1 + generators2

      frames = list(expand_xforms(g, N=depth, maxrad=maxrad))

      srcstate = state if nstates > 0 else 0
      tgtstate = -1 if nstates > 0 else 0
      makesym(frames, f'{sele} and chain A', srcstate=srcstate, tgtstate=tgtstate)

   cmd.set_view(v)
   cmd.set('suspend_updates', 'off')

def make_xtal(
   sym=None,
   sele=None,
   asym='vis and chain A',
   depth=6,
   maxrad=9e9,
   nstates=0,
   **kw,
):
   import willutil as wu

   cmd.set_state_order('MAKESYM', range(50, 0, -1))

   cmd.set('suspend_updates', 'on')
   v = cmd.get_view()

   cmd.delete('obj MAKESYM')
   cmd.delete('obj TMP_*')
   if sele is None:
      sele = cmd.get_object_list()[0]

   xtal = wu.sym.Xtal(sym)
   nchain = len(cmd.get_chains(sele))

   for istate in range(max(1, nstates)):
      state = istate + 1
      coords = cmd.get_coords(sele, state)

      coords = coords.reshape(nchain, len(coords) // nchain, 3)

      cellsize, _ = xtal.fit_coords(coords, noshift=True)
      ic(istate, cellsize)
      frames = xtal.frames(cellsize=cellsize, **kw)
      frames = [Xform(Mat(*x[:3, :3].flat), Vec(x[:3, 3])) for x in frames]
      # ic(frames)

      srcstate = state if nstates > 0 else 0
      tgtstate = -1 if nstates > 0 else 0
      makesym(frames, f'{sele} and chain A', srcstate=srcstate, tgtstate=tgtstate)

   cmd.set_view(v)
   cmd.set('suspend_updates', 'off')

def makesym(frames0, sele="all", newobj="MAKESYM", depth=None, maxrad=9e9, n=9e9, verbose=False, srcstate=0,
            tgtstate=-1,):
   v = cmd.get_view()
   if tgtstate >= 0:
      cmd.delete(newobj)
   sele = "((" + sele + ") and (not obj TMP_makesym_*))"
   selechains = cmd.get_chains(sele)
   if verbose: print(selechains)
   if not depth:
      frames = frames0
   else:
      frames = expand_xforms(frames0, N=depth, maxrad=maxrad)

   # # order on COM transform dis
   # cen = com(sele)
   # frames = sorted(frames, key=lambda x: cen.distance(x * cen))

   # make new objs
   for i, x in enumerate(frames):
      if i >= n:
         break
      # if verbose: print i, x.pretty()
      tmpname = "TMP_makesym_%i" % i
      cmd.create(tmpname, sele)
      for j, c in enumerate(selechains):
         cmd.alter('obj ' + tmpname + " and chain " + c, "chain='%s'" % ROSETTA_CHAINS[len(selechains) * i + j])
      xform(tmpname, x)
   cmd.create(newobj, "TMP_makesym_*", source_state=srcstate, target_state=tgtstate)
   cmd.delete("obj TMP_makesym_*")
   cmd.set_view(v)
   # util.cbc()

def makecx(sel='all', name="TMP", n=5, axis=Uz, chainoffset=0):
   if sel == 'all':
      for i, o in enumerate(cmd.get_object_list()):
         makecx(sel=o, name="TMP%i" % i, n=n, axis=axis)
      return
   v = cmd.get_view()
   cmd.delete("TMP__C%i_*" % n)
   chains = ROSETTA_CHAINS[chainoffset:]
   for i in range(n):
      cmd.create("TMP__C%i_%i" % (n, i), sel + " and (not TMP__C%i_*)" % n)
   for i in range(n):
      rot("obj TMP__C%i_%i" % (n, i), axis, -360.0 * float(i) / float(n))
   for i in range(n):
      cmd.alter("TMP__C%i_%i" % (n, i), "chain = '%s'" % chains[i])
   util.cbc("TMP__C*")
   # for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i),"resi=str(int(resi)+%i)"%(1000*i));
   # util.cbc("TMP__C*")
   cmd.create(name, "TMP__*")
   cmd.delete("TMP__*")
   cmd.set_view(v)
   cmd.disable(sel)
   cmd.enable(name)

def mofview():
   cmd.set('sphere_scale', 0.3)
   cmd.hide('ev')
   # cmd.show('sti')
   #cmd.show('lines', 'name n+ca+c')
   cmd.show('sti', 'resn asp+das+cys+dcs+his+dhi+glu+dgu+zn+bpy and not hydro and not name n+c+o')
   # cmd.show('sti', 'resn cys and name HG')
   cmd.show('sph', 'name ZN')

   # cmd.show('car')
   cmd.show('sti', 'name n+ca+c+cb')
   cmd.show('sph', 'name cb and not resn asp+das+cys+dcs+his+dhi+glu+dgu')

   # util.cbag('all')
   # cmd.color('green', 'name N')

   cmd.unbond('name zn', 'all')
   cmd.bond('name zn', '(not elem H+C) within 3 of name zn')

   showline(Vec(-2, -1, 1) * 20, Vec(0, 0, 0))
   showline(Vec(-1, -1, 0) * 20, Vec(0, 0, 0))
   showaxes()
   cmd.show('cgo')
   # makec3(axis=Vec(1, 1, 1))
   # cmd.hide('sti')
   # util.cbag()
   # cmd.zoom()
   # cmd.show('lin')
   # cmd.show('sti', 'resn asp+das+cys+dcs+his+dhi+glu+dgu+zn')

cmd.extend('mofview', mofview)

def diffview():
   global MOVE_UP_DOWN_FUNC

   cmd.do('remove not name n+ca+c+o; dss; util.cbc; show sph')
   # cmd.turn('y', -80)

   o = cmd.get_object_list('vis')

cmd.extend('showme', diffview)

def makedx(sel='all', n=2, name=None):
   if not name:
      name = sel.replace("+", "").replace(" ", "") + "_D%i" % n
   cmd.delete(name)
   v = cmd.get_view()
   cmd.delete("_TMP_D%i_*" % n)
   ALLCHAIN = ROSETTA_CHAINS
   chains = cmd.get_chains(sel)
   for i in range(n):
      dsel = "_TMP_D%i_%i" % (n, i)
      dsel2 = "_TMP_D%i_%i" % (n, n + i)
      cmd.create(dsel, sel + " and (not _TMP_D%i_*)" % n)
      rot(dsel, Uz, 360.0 * float(i) / float(n))
      cmd.create(dsel2, dsel)
      rot(dsel2, Ux, 180.0)
      for ic, c in enumerate(chains):
         cmd.alter("((%s) and chain %s )" % (dsel, c), "chain = '%s'" % ALLCHAIN[len(chains) * (i) + ic])
         cmd.alter("((%s) and chain %s )" % (dsel2, c), "chain = '%s'" % ALLCHAIN[len(chains) * (i + n) + ic])
   cmd.create(name, "_TMP_D*")
   util.cbc(name)
   cmd.delete("_TMP_D*")
   cmd.set_view(v)
   cmd.disable(sel)
   cmd.enable(name)

for i in range(2, 21):
   globals()['makec%i' % i] = functools.partial(makecx, n=i)
for i in range(2, 21):
   globals()['maked%i' % i] = functools.partial(makedx, n=i)

def makecxauto():
   for o in cmd.get_object_list():
      n = int(re.search("_C\\d+_", o).group(0)[2:-1])
      makecx(o, n)

def maketet(sel='all', name="TET", n=12):
   makesym(frames0=SYMTET, sele=sel, newobj=name, n=n)

def makeoct(sel='all', name="OCT", n=24):
   makesym(frames0=SYMOCT, sele=sel, newobj=name, n=n)

def makeicos(sel='all', name="ICOS", n=60):
   makesym(frames0=SYMICS, sele=sel, newobj=name, n=n)

def make_d3oct(d3, cage, cage_trimer_chain="A", depth=4, maxrad=9e9):
   if verbose:
      print(cmd.super("((" + cage + ") and (chain " + cage_trimer_chain + "))", "((" + d3 + ") and (chain A))"))
   zcagecen = com(cage + " and name ca").z
   if verbose: print(zcagecen)
   # return
   x = alignvectors(Vec(1, 1, 1), Vec(1, -1, 0), Vec(0, 0, 1), Vec(1, 0, 0))
   # if verbose: print x * Vec(1,1,1), x*Vec(1,-1,0)
   # RAD(Ux,180), RAD(Uy,120),
   G = [
      RAD(Ux, 180),
      RAD(Uz, 120),
      RAD(x * Vec(1, 0, 0), 90, Vec(0, 0, zcagecen)),
      RAD(x * Vec(1, 1, 0), 180, Vec(0, 0, zcagecen)),
   ]
   makesym(G, sele="((" + d3 + ") and ((chain A+B) and name CA))", depth=depth, maxrad=maxrad)
   cmd.show("sph", "MAKESYM")
   # cmd.disable("all")
   cmd.enable("MAKESYM")

def make_d3tet(d3, cage, cage_trimer_chain="A", depth=4, maxrad=9e9):
   if verbose:
      print(cmd.super("((" + cage + ") and (chain " + cage_trimer_chain + "))", "((" + d3 + ") and (chain A))"))
   zcagecen = com(cage + " and name ca").z
   if verbose: print(zcagecen)
   # return
   x = alignvectors(Vec(1, 1, 1), Vec(1, -1, 0), Vec(0, 0, 1), Vec(1, 0, 0))
   # if verbose: print x * Vec(1,1,1), x*Vec(1,-1,0)
   # RAD(Ux,180), RAD(Uy,120),
   G = [
      RAD(Ux, 180),
      RAD(Uz, 120),
      RAD(x * Vec(1, 0, 0), 180, Vec(0, 0, zcagecen)),
   ]
   makesym(G, sele="((" + d3 + ") and ((chain A+B) and name CA))", depth=depth, maxrad=maxrad)
   cmd.show("sph", "MAKESYM")
   # cmd.disable("all")
   cmd.enable("MAKESYM")

def symfit(sym='icos', sele='visible'):
   import willutil as wu
   import numpy as np
   chains = pymol.cmd.get_chains(sele)
   coords = np.stack(cmd.get_coords(f'({sele} and name CA and chain {c})') for c in chains)
   assert len(coords) == 4

   # ic(coords[0].shape)
   # ic(coords[1].shape)
   rms1, _, xfit1 = wu.hrmsfit(coords[1], coords[0])
   rms2, _, xfit2 = wu.hrmsfit(coords[2], coords[0])
   rms3, _, xfit3 = wu.hrmsfit(coords[3], coords[0])
   # ic(rms1,rms2,rms3)

   frames = np.stack([np.eye(4), xfit1, xfit2, xfit3])
   symfit = wu.compute_symfit('icos', frames)
   # ic(symfit.xfit)

   crd = cmd.get_coords(sele)
   crd = wu.hxform(symfit.xfit, crd)
   crd = wu.unhomog(crd)
   cmd.load_coords(crd, sele)

cmd.extend('symfit', symfit)
