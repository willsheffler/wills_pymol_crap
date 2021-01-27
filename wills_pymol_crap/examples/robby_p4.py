import pymol
from pymol import cmd
from wills_pymol_crap import symgen
from copy import deepcopy
from wills_pymol_crap.xyzMath import Vec
from wills_pymol_crap.symgen_test import test_xtal
from wills_pymol_crap.symgen import (
   SymElem,
   generate_sym_trie,
   BuildCGO,
   ComponentCenterVisitor,
   CountFrames,
)

def main():
   cellsize = 100

   p32gen = [
      SymElem("C4", axis=Vec(0, 0, 1), cen=Vec(0.0, 0.0, 0.0) * cellsize),
      SymElem("C2", axis=Vec(0, 0, 1), cen=Vec(0.5, 0.0, 0.0) * cellsize),
   ]
   robby_make_symfiles(generators=p32gen, tag='P42_C4_C2', depth=6, cellsize=cellsize)

def robby_make_symfiles(generators, tag, depth, cellsize):
   cmd.delete(f'all')

   print(f'making symdef files for {tag:}')

   test_xtal(
      generators,
      tag=tag,
      cell=cellsize,
      depth=depth,
      mindepth=2,  # does all from [mindepth..depth]
      symdef=True,
      one_component=True,
      symdef_scale=1 / cellsize / 1000,
      generic_names=True,
   )
   cmd.zoom()

   fname = f'{tag}.pse'
   print(f'saving session to {fname:}')
   cmd.save(fname)

   cmd.delete(all)

def robby_make_P4212_C4_C2(cell=80, **kw):
   pass

if __name__ == '__main__':
   main()