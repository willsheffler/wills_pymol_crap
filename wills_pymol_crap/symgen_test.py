import random

from symgen_classes import (
   SymElem,
   generate_sym_trie,
   BuildCGO,
   ComponentCenterVisitor,
   CountFrames,
)
from pymol import cmd

def test_xtal(
   G,
   *,
   cell,
   depth=4,
   mindepth=0,
   symdef=1,
   shownodes=1,
   verbose=False,
   showcell=False,
   **kw,
):
   print('test_xtal')
   v = cmd.get_view()
   CEN = [g.cen for g in G]
   FN = list()
   tag = "test" if not "tag" in kw else kw["tag"]
   for d in range(mindepth, depth + 1):
      symtrie = generate_sym_trie(G, depth=d)
      # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
      nodes = []
      if "component_pos" in list(kw.keys()):
         raise NotImplementedError("component_pos is no longer used")
         # nodes = kw["component_pos"][:1]
      buildcgo = BuildCGO(nodes=nodes, label=tag + "_DEPTH%i" % d, cell=cell, **kw)
      symtrie.visit(buildcgo)
      buildcgo.show()
      if shownodes:
         find_nodes = ComponentCenterVisitor(symelems=G, label=tag + "_NODES%i" % d, **kw)
         symtrie.visit(find_nodes)
         FN.append(find_nodes)
         if symdef:
            sdef_string = FN[-1].make_symdef(**kw)
            if verbose:
               print("==================== SYMDEF (dump to " + tag + "_" + str(d) + ".sym) ====================")
               print(sdef_string)
               print("=====================================================================")
            with open(tag + "_" + str(d) + ".sym", "w") as out:
               out.write(sdef_string)
            if "symdef_check" in kw and kw["symdef_check"]:
               sdef = RosettaSymDef()
               sdef.parse(sdef_string)
               sdef.show("SYMDEF_" + tag)
   for fn in FN:
      fn.show(**kw)  # dumb order hack for pymol up/dn
   cmd.disable("all")
   cmd.enable(tag + "_DEPTH%i" % (depth))
   cmd.enable(tag + "_NODES%i" % (depth))
   count = CountFrames()
   symtrie.visit(count)
   if verbose: print("N Frames:", count.count)

   if showcell:
      cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))

   cmd.set_view(v)

# I213
#    AXS = [Vec(1, 1, 1), Vec(1, 0, 0)]
#    CEN = [cell * Vec(0, 0, 0), cell * Vec(0, 0, 0.25)]

def test_I4132(cell=100, tag='I4132', **kw):
   """
delete all; run /home/sheffler/pymol3/symgen_classes.py; run /home/sheffler/pymol3/symgen.py; run /home/sheffler/pymol3/symgen_test.py; test_I4132(depth=7,shownodes=0,cell=200,maxrad=180); run /home/sheffler/pymol3/misc/G222.py; gyroid(200,r=180)
   """
   # delete all; run /Users/sheffler/pymol/symgen.py; test_I4132( cell=150,
   # depth=7, shownodes=0, maxrad=150 );  gyroid(150,r=150,c=Vec(75,75,75))
   G = [
      # SymElem("C3", cen=cell * Vec(0.5, 0.5, 0.5), axis=Vec(1, 1, 1), col=(1, 0.5, 0.5)),
      # SymElem("C2", cen=cell * Vec(1.125 / 2, 0.500, 0.625), axis=Vec(1, 0, 0), col=(0.5, 1, 0.5)),
      # SymElem("C2", cen=cell * Vec(0.875 / 2, 0.500, 0.375), axis=Vec(1, 0, 0), col=(0.5, 0.5, 1)),
      # SymElem("C3", axis=Vec(1, 1, 1)),
      # SymElem("C2", axis=Vec(1, 1, 0), cen=cell * Vec(-1, 1, 1) / 8.0),
      # SymElem("C2", axis=Vec(1, 1, 0), cen=cell * Vec(1, -1, -1) / 8.0),
      # SymElem( "C2", axis=Vec(1,1,0), cen=cell*Vec( 1,-1,-1)/8.0, col=(1,1,0) ),
      SymElem("D3", cen=cell * Vec(0.125, 0.125, 0.125), axis=Vec(1, 1, 1), axis2=Vec(1, -1, 0), col=(0, 1, 0)),
      # SymElem("D2", cen=cell * Vec(0.125, 0.000, 0.250), axis=Vec(1, 0, 0), axis2=Vec(0, -1, 1), col=(0, 1, 1)),
      # SymElem("D3", cen=cell * Vec(-0.125, -0.125, -0.125), axis=Vec(1, 1, 1), axis2=Vec(1, -1, 0), col=(1, 0, 0)),
      SymElem("D2", cen=cell * Vec(-0.125, 0.000, -0.250), axis=Vec(1, 0, 0), axis2=Vec(0, -1, 1), col=(1, 1, 0)),
   ]
   test_xtal(G, cell=cell, tag=tag, origin=cell * Vec(0.5, 0.5, 0.5), showshape=0, **kw)
   # cube( cell*Vec(-0.5,-0.5,-0.5), cell*Vec(0.5,0.5,0.5) )
   # cube(cell * Vec(-0, -0, -0), cell * Vec(1, 1, 1))

def test_I213(cell=70, **kw):
   """
delete all; run ~/pymol3/symgen.py; run ~/pymol3/symgen_test.py; test_I213(depth=6,mindepth=6,cell=100)
   """
   # AXS = [Vec(1, 1, 1), Vec(1, 1, -1)]
   # CEN = [cell * Vec(0, 0, 0), cell * Vec(0.5, 0, 0.0)]
   G = [
      SymElem("C2", axis=Vec(1, 0, 0), cen=Vec(0, 0, 0.25) * cell, col=[1, 0.7, 0.0]),
      SymElem("C3", axis=Vec(1, 1, 1) * 0.57735, cen=Vec(0, 0, 0) * cell, col=[0.1, 0.5, 1]),

      # SymElem("C2", axis=Vec(0, 1, 0), cen=Vec(-1, 0, -2)/4 * cell, col=[1, 0, 0]),
      # SymElem("C3", axis=Vec(1, -1, 1), cen=Vec(2,1,-1)/6 * cell, col=[0, 1, 0]),
   ]

   test_xtal(
      G,
      cell,
      tag="I213",
      origin=cell * Vec(0.0, 0.0, 0.0),
      showshape=0,
      symdef=0,
      make_symdef=0,
      showlinks=False,
      shownodes=False,
      radius=0.5,
      bbox=[Vec(cell, cell, cell) * -0.1, Vec(cell, cell, cell) * 1.1],
      # length=40,
      **kw,
   )
   showcube(Vec(0, 0, 0), Vec(1, 1, 1) * cell)

def test_P213(cell=70, **kw):
   # P213   C3 -0.57735,0.57735,0.57735   -  -0.5,-0.5,0 C3 0.57735,0.57735,0.57735 -  -0.333333,-0.333333,0.666667  70.5288  0.353553 -
   'delete all; run ~/pymol3/pymol_util.py; ~/pymol3/symgen.py; run ~/pymol3/symgen_test.py; test_P213(depth=2)'
   # AXS = [Vec(1, 1, 1), Vec(1, 1, -1)]
   # CEN = [cell * Vec(0, 0, 0), cell * Vec(0.5, 0, 0.0)]
   AXS, CEN = [], []
   AXS.append(Vec(-0.57735, 0.57735, 0.57735))
   CEN.append(Vec(-0.5, -0.5, 0) * cell)
   AXS.append(Vec(0.57735, 0.57735, 0.57735))
   # CEN.append(Vec(0.333333, -0.333333, 0.666667) * cell) # straight from franks table... WTF???
   CEN.append(Vec(0, 0, 0) * cell)
   G = [
      SymElem("C3", axis=AXS[0], cen=CEN[0], col=[1, 0.7, 0.0]),
      SymElem("C3", axis=AXS[1], cen=CEN[1], col=[0.1, 0.5, 1]),
   ]
   test_xtal(
      G,
      cell=cell,
      tag="P213",
      origin=cell * Vec(0.0, 0.0, 0.0),
      showshape=0,
      symdef=1,
      make_symdef=1,
      showlinks=1,
      radius=0.5,
      # length=40,
      **kw,
   )

def test_P4132(cell=100, **kw):
   # P4132 C2 0,-0.707107,0.707107 -  0.375,-0.125,-0.125  C3 0.57735,-0.57735,0.57735   -  0.333333,0.166667,-0.166667   35.2644  0.176777 -
   'delete all; run ~/pymol3/symgen.py; run ~/pymol3/symgen_test.py; test_P4132(depth=2)'
   AXS, CEN = [], []
   AXS.append(Vec(0, -0.707107, 0.707107))
   CEN.append(Vec(0.375, -0.125, -0.125) * cell)
   # AXS.append(Vec(0.57735, -0.57735, 0.57735)),
   # CEN.append(Vec(0.333333, 0.166667, -0.166667) * cell)
   AXS.append(Vec(0.57735, -0.57735, 0.57735)),
   CEN.append(Vec(0.333333, 0.166667, -0.166667) * cell)
   # CEN.append(Vec(0, 0, 0) * cell)
   G = [
      SymElem("C2", axis=AXS[0], cen=CEN[0], col=[1, 0.7, 0.0]),
      SymElem("C3", axis=AXS[1], cen=CEN[1], col=[0.1, 0.5, 1]),
   ]
   test_xtal(
      G,
      cell,
      tag="P4132",
      origin=cell * Vec(0.0, 0.0, 0.0),
      showshape=0,
      symdef=0,
      showlinks=0,
      radius=0.5,
      # length=40,
      **kw,
   )
   cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))

def test_P4132_old(depth=8, cell=50, maxrad=80):
   # delete all; run /home/sheffler/pymol3/symgen_test.py; test_P4132(depth=8, cell=50, maxrad=80)
   # **** P4132 ****
   # C2 and C3 at angle = 35.2644 offset = 0.176777
   #     C2 axis=[-0.707107,0.707107,0]  origin=[-0.125,-0.125,-0.125]
   #     C3 axis=[0.57735,-0.57735,0.57735]  origin=[0,-0.5,-0.5]
   AXS = [Vec(-1, 1, 0), Vec(1, -1, 1)]
   CEN = [cell * Vec(1, 1, 1) / -8.0, cell * Vec(0, 1, 1) / -2.0]
   # AXS = [ Vec(-1, 1, 0) ,
   #         Vec( 1, 1, 1) ]
   # CEN = [ cell * Vec(1,1,1)/-8.0,
   #         cell * Vec(0,0,0) ]
   G = [SymElem("C2", axis=AXS[0], cen=CEN[0]), SymElem("C3", axis=AXS[1], cen=CEN[1])]
   symtrie = generate_sym_trie(G, depth=depth)
   # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
   cencell = cell / 2.0 * Vec(1, 1, 1)
   buildcgo = BuildCGO(
      nodes=[CEN[1] + randnorm() * 5.0, CEN[0] + randnorm() * 8.0],
      origin=cencell,
      maxrad=maxrad,
      showlinks=False,
      showelems=True,
   )
   symtrie.visit(buildcgo)
   buildcgo.show()
   cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
   for g in G:
      print("show", g)
      g.show(radius=2.0, sphereradius=4.0)
   return AXS, CEN

# def test_P23_D2TET(depth=6, cell=60, **kw):
#    G = [
#       # SymElem("C2",axis=Vec( 1,0,0),cen=Vec(    0    ,cell/2.0, cell/4.0)),
#       # SymElem("C3",axis=Vec(-1,1,1),cen=Vec(-cell/6.0,cell/6.0,-cell/3.0)),
#       # SymElem("C2",axis=Vec(1,1,0),cen=cell*Vec(0,0.0,0.0)),
#       # SymElem("C3",axis=Vec(-1,1,1),cen=Vec(0,0,0)),
#       # SymElem("C2",axis=Vec(0,1,0)),
#       SymElem("T", axis=Vec(1, 1, 1)),
#       # SymElem("C2",axis=Vec(0,1,0),cen=Vec(1,1,1)*cell/2.0), # other T
#       # SymElem("C3",axis=Vec(1,1,1),cen=Vec(1,1,1)*cell/2.0), # other T
#       # SymElem("C4",axis=Vec(1,0,0),cen=Vec(1,1,1)*cell/2.0),
#       SymElem("D2", cen=Vec(0, 0, cell / 2.0)),
#       # SymElem("D2",cen=Vec(0,cell/2.0,cell/2.0)), # other D2
#    ]
#    # for elem in G: print elem
#    symtrie = generate_sym_trie(G, depth=depth)
#    buildcgo = BuildCGO(nodes=[], origin=Vec(0.5, 0.5, 0.5) * cell, **kw)
#    # buildcgo = BuildCGO( nodes=[ Vec(6,6,30), ] )
#    symtrie.visit(buildcgo)
#    buildcgo.show()
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
#    for g in G:
#       print("show", g)
#       g.show(radius=2.0, sphereradius=4.0)

# def test_F23_TETC2(cell=100, **kw):
#    #   )
#    G = [
#       SymElem("T", axis=Vec(1, 1, 1)),
#       SymElem("C2", axis=Vec(0, 0, 1), cen=Vec(cell / 4.0, cell / 4.0, 0)),
#    ]
#    test_xtal(G, cell, tag="F23_TC2", **kw)

# def test_F23_TETC3(cell=100, **kw):
#    G = [
#       SymElem("T", axis=Vec(1, 1, 1)),
#       SymElem("C3", axis=Vec(-1, 1, 1), cen=Vec(cell / 3, cell / 6, cell / 6)),
#    ]
#    test_xtal(G, cell, tag="F23_TC3", **kw)
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))

# def test_F23_TT(cell=100, **kw):
#    G = [SymElem("T"), SymElem("T", cen=Vec(cell / 2, cell / 2, 0))]
#    test_xtal(G, cell, tag="F23_TT", **kw)
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))

# def test_D4OCT(depth=4):
#    CELL = 30
#    G = [
#       # SymElem("C2",axis=Vec( 1,0,0),cen=Vec(    0    ,CELL/2.0, CELL/4.0)),
#       # SymElem("C3",axis=Vec(-1,1,1),cen=Vec(-CELL/6.0,CELL/6.0,-CELL/3.0)),
#       # SymElem("C2",axis=Vec(1,1,0),cen=CELL*Vec(0,0.0,0.0)),
#       # SymElem("C3",axis=Vec(-1,1,1),cen=Vec(0,0,0)),
#       SymElem("C2", axis=Vec(1, 1, 0), cen=Vec(0, 0, CELL)),
#       SymElem("C3", axis=Vec(1, 1, 1), cen=Vec(0, 0, CELL)),
#       SymElem("C4", axis=Vec(1, 0, 0), cen=Vec(0, 0, CELL)),
#       SymElem("C4", axis=Vec(0, 0, 1)),
#       SymElem("C2", axis=Vec(1, 0, 0)),
#       SymElem("C2", axis=Vec(0, 1, 0)),
#    ]
#    # for elem in G: print elem
#    symtrie = generate_sym_trie(G, depth=depth)
#    # buildcgo = BuildCGO( nodes=[ Vec(2,4,3), Vec(-6,-2,35), ] )
#    buildcgo = BuildCGO(nodes=[Vec(-6, -2, CELL + 5), Vec(2, 4, 3)])
#    # buildcgo = BuildCGO( nodes=[ Vec(-6,-2,35), ] )
#    # buildcgo = BuildCGO( nodes=[ Vec(2,4,3), ] )
#    symtrie.visit(buildcgo)
#    buildcgo.show()

def test_I432_OD3(cell=120, **kw):
   X = alignvectors(Vec(1, 1, 1), Vec(1, 0, -1), Vec(0, 0, 1), Vec(1, 0, 0))
   G = [
      SymElem("D3", cen=cell * Vec(0.0, 0.0, 0.0)),
      SymElem("O", cen=cell * Vec(0.0, 0.0,
                                  sqrt(3.0) / 4.0), input_xform=X),
   ]
   cube(cell * Vec(-0.25, -0.25, -0.25), cell * Vec(0.75, 0.75, 0.75), xform=X)
   test_xtal(G, cell, tag="I432_OD3", one_component=True, **kw)

def test_octa(cell=100, **kw):
   X = alignvectors(Vec(1, 1, 1), Vec(1, 0, -1), Vec(0, 0, 1), Vec(1, 0, 0))
   G = [
      SymElem("D3", axis=Vec(1, 0, 0), axis2=Vec(0, -1, 1), cen=cell * Vec(0.0, 0.25, 0.25)),
      SymElem("O", cen=cell * Vec(0.0, 0.0, 0), input_xform=X),
   ]
   cube(cell * Vec(-0.25, -0.25, -0.25), cell * Vec(0.75, 0.75, 0.75))
   test_xtal(G, cell, tag="I432_OD3", one_component=True, **kw)

# F432  O  0.57735, 0.57735,0.57735   0,0,0       D2  1,0,0 0,-0.707107,0.707107 0,0.25,0.25
# F432  O  0.57735,-0.57735,0.57735   0,-0.5,-0.5 D2 -0.707107,0.707107,0 0.707107,0.707107,0  -1,0,0
# F432  O  0.57735,-0.57735,0.57735   0,-0.5,-0.5 D3 -0.707107,0.707107,0 0.57735,0.57735,0.57735
# I432  O  0.57735, 0.57735,0.57735   0,0,0       D3  0,-0.707107,0.707107 0.57735,0.57735,0.57735
# F432  O  0.57735,-0.57735,0.57735   0,-0.5,-0.5 D4  0.707107,0.707107,0  0,0,1 0,0,-0.5 90 0.204124 -

# def test_P23_TT(cell=100, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P23_TD2B( depth=2, cell=200,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("C3", axis=Vec(1, -1, -1), cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C3", axis=Vec(1, -1, 1), cen=cell * Vec(0.25, 0.25, 0.0)),
#       # SymElem( "T" , cen=cell*Vec(0.25,0.25,0.0) ),
#    ]
#    test_xtal(G, cell, tag="P23_TT", **kw)
#    cube(cell * Vec(1, 1, 1) * -0.5, cell * Vec(1, 1, 1) * 0.5)

# def test_P432_OD4(cell=80, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("O", cen=cell * Vec(0.0, 0.0, 0.5)),
#    ]
#    test_xtal(G, cell, tag="P432_OD4", **kw)

# def test_P4212_C4_C2(cell=80, **kw):
#    # delete all; run /home/sheffler/pymol3/symgen_test.py; test_P4212_C4_C2( depth=6, cell=1 )

#    generators = [
#       SymElem("C4", axis=Vec(0, 0, 1), cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C2", axis=Vec(0, 1, 0), cen=cell * Vec(0.5, 0.0, 0.0)),
#    ]
#    test_xtal(generators, cell, tag="P4212_C4_C2", **kw)

# def test_P4m_44(cell=80, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("D4", cen=cell * Vec(0.5, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P432_OD4", **kw)

# def test_P23_TD2A(cell=200, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P23_TD2B( depth=2, cell=200,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("T", cen=cell * Vec(0.0, 0.0, 0.5)),
#    ]
#    component_pos = [Vec(-8, -7, 6), Vec(-7, 0, -11)]
#    test_xtal(G, cell, component_pos=component_pos, tag="P23_TD2A", **kw)
#    cube(cell * Vec(0, 0, -0.5), cell * Vec(1, 1, 0.5))

# def test_P23_TC3(cell=200, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P23_TD2B( depth=2, cell=200,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("C3", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("T", cen=cell * Vec(0.5, 0.5, 0.5)),
#    ]
#    test_xtal(G, cell, tag="P23_TC3", **kw)
#    cube(cell * Vec(0, 0, 0), cell * Vec(1, 1, 1))

# # run /Users/sheffler/pymol/symgen.py; delete all; test_P23_TD2B( depth=4 )

# def test_P23_TD2B(depth=3, cell=150, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P23_TD2B( depth=2, cell=200,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("T", cen=cell * Vec(0.0, 0.5, 0.5)),
#    ]
#    test_xtal(G, cell, tag="P23_TD2B", **kw)
#    cube(cell * Vec(0, -0.5, -0.5), cell * Vec(1, 0.5, 0.5))

# def test_F432_TD2(cell=150, **kw):
#    # delete all; run ~/pymol/symgen.py; test_F432_TD2( depth=2, cell=200,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),  # , axis2=Vec(1,1,0) ),
#       SymElem("T", cen=cell * Vec(0.0, 0.0, 0.25), input_xform=RAD(Uz, 45)),
#       # SymElem( "O" , cen=cell*Vec(0.5,0.5,0.0) ),
#    ]
#    component_pos = [Vec(2, 4, 8), Vec(2, 5, -12), Vec(4, 6, 8)]
#    test_xtal(G, cell, component_pos=component_pos, tag="F432_TD2", **kw)
#    cube(cell * Vec(-0.5, -0.5, -0.25), cell * Vec(0.5, 0.5, 0.75), xform=RAD(Uz, 45))

# def test_F432_OTD2(depth=3, cell=80, **kw):
#    # delete all; run /Users/sheffler/pymol/symgen.py; test_P23(depth=1,mindepth=1)
#    # ODT 96
#    # TDO 96
#    G = [
#       SymElem("O", cen=cell * Vec(0.5, 0.5, 0.0)),
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0), axis2=Vec(1, 1, 0)),
#       SymElem("T", cen=cell * Vec(0.0, 0.0, -0.5)),
#    ]
#    component_pos = [Vec(2, 4, 8), Vec(2, 5, -12), Vec(4, 6, 8)]
#    extranodes = [cell * Vec(0.0, 0.5, 0.0)]
#    test_xtal(G, cell, component_pos=component_pos, tag="F432_TD2", extranodes=extranodes, **kw)
#    cube(cell * Vec(-1, -1, -0.5), cell * Vec(1, 1, 1.5))

# def test_I432_OD4(cell=80, **kw):
#    # delete all; run /Users/sheffler/pymol/symgen.py; test_F432_OD4( depth=3, cell=200, symdef_scale=0.000001, generic_names=1 )
#    # G = [ SymElem( "D2", cen=cell*Vec(0.0,0.0,0.0), axis2=Vec(1,1,0) ),
#    #       # SymElem( "T" , cen=cell*Vec(0.0,0.0,-0.5) ),
#    #       SymElem( "O" , cen=cell*Vec(0.5,0.5,0.0) ),
#    #     ]
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       # SymElem( "T" , cen=cell*Vec(0.0,0.0,-0.5) ),
#       SymElem("O", cen=cell * Vec(0.0, 0.0, 0.5), input_xform=RAD(Uz, 45)),
#    ]
#    component_pos = [Vec(2, 4, 8), Vec(2, 5, -12), Vec(4, 6, 8)]
#    test_xtal(G, cell, component_pos=component_pos, tag="I432_OD4", **kw)
#    cube(cell * Vec(-1, -1, -0.5), cell * Vec(1, 1, 1.5))

# def test_P442(cell=80, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P4( depth=4, cell=100,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("C4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C2", cen=cell * Vec(0.5, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="test_P4_42", **kw)

# def test_P444(cell=100, **kw):
#    G = [
#       SymElem("C4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       # moved lattice WRP C2!!!!!
#       SymElem("C4", cen=cell * Vec(0.5, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="test_P4_44", **kw)

# def test_P4m44(cell=100, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       # moved lattice WRP C2!!!!!
#       SymElem("D4", cen=cell * Vec(0.5, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="test_P4_44", **kw)

# def test_I4132_D(cell=100, **kw):
#    # delete all; test_I4132(depth=7,shownodes=0,cell=200,maxrad=180); run /Users/sheffler/pymol3/misc/G222.py; gyroid(200,r=180)
#    # delete all; run /Users/sheffler/pymol/symgen.py; test_I4132( cell=150,
#    # depth=7, shownodes=0, maxrad=150 );  gyroid(150,r=150,c=Vec(75,75,75))
#    G = [
#       SymElem(
#          "D3",
#          cen=cell * Vec(0.625, 0.625, 0.625),
#          axis=Vec(1, 1, 1),
#          axis2=Vec(1, -1, 0),
#          col=(0, 1, 0),
#       ),
#       SymElem(
#          "D2",
#          cen=cell * Vec(0.625, 0.500, 0.750),
#          axis=Vec(1, 0, 0),
#          axis2=Vec(0, -1, 1),
#          col=(0, 1, 1),
#       ),
#       SymElem(
#          "D3",
#          cen=cell * Vec(0.375, 0.375, 0.375),
#          axis=Vec(1, 1, 1),
#          axis2=Vec(1, -1, 0),
#          col=(1, 0, 0),
#       ),
#       SymElem(
#          "D2",
#          cen=cell * Vec(0.375, 0.500, 0.250),
#          axis=Vec(1, 0, 0),
#          axis2=Vec(0, -1, 1),
#          col=(1, 1, 0),
#       ),
#       #       SymElem( "C3", axis=Vec(1,1,1) ),
#       # SymElem( "C2", axis=Vec(1,1,0), cen=cell*Vec(-1, 1, 1)/8.0 ),
#       #       SymElem( "C2", axis=Vec(1,1,0), cen=cell*Vec( 1,-1,-1)/8.0 ),
#       # SymElem( "C2", axis=Vec(1,1,0), cen=cell*Vec( 1,-1,-1)/8.0, col=(1,1,0) ),
#       # SymElem( "D3", cen=cell*Vec(0.125,0.125,0.125), axis=Vec(1,1,1), axis2=Vec(1,-1,0), col=(0,1,0) ),
#       # SymElem( "D2", cen=cell*Vec(0.125,0.000,0.250), axis=Vec(1,0,0), axis2=Vec(0,-1,1), col=(0,1,1) ),
#       # SymElem( "D3", cen=cell*Vec(-0.125,-0.125,-0.125), axis=Vec(1,1,1), axis2=Vec(1,-1,0), col=(1,0,0) ),
#       # SymElem( "D2", cen=cell*Vec(-0.125,0.000,-0.250), axis=Vec(1,0,0), axis2=Vec(0,-1,1), col=(1,1,0) ),
#    ]
#    test_xtal(G, cell, tag="I4132", origin=cell * Vec(0.5, 0.5, 0.5), showshape=0, **kw)
#    # cube( cell*Vec(-0.5,-0.5,-0.5), cell*Vec(0.5,0.5,0.5) )
#    cube(cell * Vec(-0, -0, -0), cell * Vec(1, 1, 1))

# cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))

# def test_T32(cell=100, **kw):
#    print(
#       "========================================== T32 ==========================================")
#    G = [
#       SymElem("C3", axis=Vec(1, 1, 1)),
#       SymElem("C2", axis=Vec(1, 0, 0), cen=Vec(0.00001, 0, 0)),
#    ]
#    component_pos = [Vec(8, 6, 5), Vec(7, 2, 1)]
#    test_xtal(G, cell, component_pos=component_pos, tag="T32", origin=cell * Vec(0.0, 0.0, 0.0),
#              showshape=0, **kw)

# def test_I213(depth=16, cell=100, maxrad=9e9):
#    # C3 and C2 at angle = 54.7356 offset = 0.176777
#    #     C3 axis=[-0.57735,0.57735,0.57735]  origin=[-0.166667,0.166667,-0.333333]
#    #     C2 axis=[1,0,0]  origin=[0,0.5,0.25]
#    # AXS = [ Vec(-1,1,1),
#    #         Vec( 1,0,0) ]
#    # CEN = [ cell * Vec(-1,1,-2)/6.0 ,
#    #         cell * Vec(0,0.5,0.25) ]
#    AXS = [Vec(1, 1, 1), Vec(1, 0, 0)]
#    CEN = [cell * Vec(0, 0, 0), cell * Vec(0, 0, 0.25)]
#    G = [SymElem("C3", axis=AXS[0], cen=CEN[0]), SymElem("C2", axis=AXS[1], cen=CEN[1])]
#    symtrie = generate_sym_trie(G, depth=depth)
#    # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
#    nodes = [
#       CEN[0] + projperp(AXS[0], randnorm()) * 8.0,
#       CEN[1] + projperp(AXS[1], randnorm()) * 8.0,
#    ]
#    buildcgo = BuildCGO(nodes=nodes, maxrad=maxrad, origin=cell * Vec(0.5, 0.5, 0.5),
#                        showlinks=True)
#    symtrie.visit(buildcgo)
#    buildcgo.show()
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
#    for g in G:
#       print("show", g)
#       g.show(radius=2.0, sphereradius=4.0)
#    return AXS, CEN

def test_P4132_old(depth=8, cell=50, maxrad=80):
   # delete all; run /home/sheffler/pymol3/symgen_test.py; test_P4132(depth=8, cell=50, maxrad=80)
   # **** P4132 ****
   # C2 and C3 at angle = 35.2644 offset = 0.176777
   #     C2 axis=[-0.707107,0.707107,0]  origin=[-0.125,-0.125,-0.125]
   #     C3 axis=[0.57735,-0.57735,0.57735]  origin=[0,-0.5,-0.5]
   AXS = [Vec(-1, 1, 0), Vec(1, -1, 1)]
   CEN = [cell * Vec(1, 1, 1) / -8.0, cell * Vec(0, 1, 1) / -2.0]
   # AXS = [ Vec(-1, 1, 0) ,
   #         Vec( 1, 1, 1) ]
   # CEN = [ cell * Vec(1,1,1)/-8.0,
   #         cell * Vec(0,0,0) ]
   G = [SymElem("C2", axis=AXS[0], cen=CEN[0]), SymElem("C3", axis=AXS[1], cen=CEN[1])]
   symtrie = generate_sym_trie(G, depth=depth)
   # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
   cencell = cell / 2.0 * Vec(1, 1, 1)
   buildcgo = BuildCGO(
      nodes=[CEN[1] + randnorm() * 5.0, CEN[0] + randnorm() * 8.0],
      origin=cencell,
      maxrad=maxrad,
      showlinks=False,
      showelems=True,
   )
   symtrie.visit(buildcgo)
   buildcgo.show()
   cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
   for g in G:
      print("show", g)
      g.show(radius=2.0, sphereradius=4.0)
   return AXS, CEN

# def test_F432(depth=6, cell=100, maxrad=90):
#    # C3 and D2 at angle = 35.2644 offset = 0
#    #     C3 axis=[0.57735,0.57735,0.57735]  origin=[0,0,0]
#    #     D2 axis=[1,0,0]  axis2=[0,-0.707107,0.707107]  origin=[0,0.25,0.25]
#    AXS = [Vec(1, 1, 1), Vec(1, 0, 0)]
#    CEN = [cell * Vec(0, 0, 0), cell * Vec(0, 0.25, 0.25)]
#    G = [
#       SymElem("C2", axis=Vec(1, 1, 0), cen=Vec(0, 0, 0)),
#       SymElem("C3", axis=AXS[0], cen=CEN[0]),
#       SymElem("C4", axis=Vec(1, 0, 0), cen=Vec(0, 0, 0)),
#       SymElem("D2", axis=AXS[1], axis2=Vec(0, -1, 1), cen=CEN[1]),
#    ]
#    symtrie = generate_sym_trie(G, depth=depth)
#    # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
#    nodes = [
#       CEN[0] + projperp(AXS[0], randnorm()) * 8.0,
#       CEN[1] + projperp(AXS[1], randnorm()) * 8.0,
#    ]
#    buildcgo = BuildCGO(nodes=nodes, maxrad=maxrad, origin=cell * Vec(0.5, 0.5, 0.5),
#                        showlinks=False)
#    symtrie.visit(buildcgo)
#    buildcgo.show()
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
#    for g in G:
#       print("show", g)
#       g.show(radius=2.0, sphereradius=4.0)
#    return AXS, CEN

# def test_quasi(depth=8, cell=20.0, maxrad=9e9):
#    G = [
#       SymElem("C2", axis=Vec(1, 0, 0)),
#       SymElem("C2", axis=Vec(0, 1, 0)),
#       SymElem("C2", axis=Vec(0, 0, 1)),
#       SymElem("C3", axis=Vec(+1, +1, +1)),
#       SymElem("C3", axis=Vec(+1, -1, -1)),
#       SymElem("C3", axis=Vec(-1, +1, -1)),
#       SymElem("C3", axis=Vec(-1, -1, +1)),
#       SymElem("C2", axis=Vec(2, -1, -1), cen=cell * Vec(1, 1, 1)),
#    ]
#    nodes = []
#    symtrie = generate_sym_trie(G, depth=depth)
#    symtrie.visit(print_node)
#    print(symtrie)
#    buildcgo = BuildCGO(
#       nodes=nodes,
#       maxrad=maxrad,
#       origin=cell * Vec(0.5, 0.5, 0.5),
#       showlinks=False,
#       showelems=True,
#    )
#    symtrie.visit(buildcgo)
#    buildcgo.show()

# def test_P432_D4_C3(cell=80, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.5, 0.5, 0.0)),
#       SymElem("C3", cen=cell * Vec(0.25, 0.25, 0.25), axis=Vec(1, 1, 1)),
#    ]
#    test_xtal(G, cell, tag="P432_D4_C3", **kw)
#    cube(cell * Vec(-0.0, -0.0, -0), cell * Vec(1, 1, 1))

# def test_I432(cell=50, **kw):
#    # **** I432 ****
#    # C2 and D3 at angle = 35.2644 offset = 0.353553
#    #     C2 axis=[-0.707107,0.707107,0]  origin=[0.25,0.25,0.25]
#    # D3 axis=[-0.57735,0.57735,0.57735]  axis2=[0.707107,0.707107,0]
#    # origin=[0,0,0]
#    G = [
#       SymElem("D3", axis=Vec(-1, 1, 1), axis2=Vec(1, 1, 0)),
#       SymElem("C2", axis=Vec(-1, 1, 0), cen=Vec(1, 1, 1) / 4.0 * cell),
#    ]
#    test_xtal(G, cell, tag="I43", origin=cell * Vec(0.5, 0.5, 0.5), showshape=0, **kw)
#    # cube( cell*Vec(-0.5,-0.5,-0.5), cell*Vec(0.5,0.5,0.5) )
#    cube(cell * Vec(-0, -0, -0), cell * Vec(1, 1, 1))

#    # symtrie = generate_sym_trie(G, depth=depth)
#    # # buildcgo = BuildCGO( nodes=[ CEN1+Vec(2,3,4), CEN2+Vec(2,4,3), ] )
#    # nodes = []
#    # buildcgo = BuildCGO(
#    #     nodes=nodes, origin=cell * Vec(0.5, 0.5, 0.5), showlinks=False, **kw
#    # )
#    # symtrie.visit(buildcgo)
#    # buildcgo.show()
#    # cube(Vec(0, 0, 0), cell * Vec(1, 1, 1))
#    # for g in G:
#    #     print("show", g)
#    #     g.show(radius=2.0, sphereradius=4.0)

# def test_I4132_half(cell=100, **kw):
#    # delete all; test_I4132(depth=7,shownodes=0,cell=200,maxrad=180); run /Users/sheffler/pymol3/misc/G222.py; gyroid(200,r=180)
#    # delete all; run /Users/sheffler/pymol/symgen.py; test_I4132( cell=150,
#    # depth=7, shownodes=0, maxrad=150 );  gyroid(150,r=150,c=Vec(75,75,75))
#    G = [
#       SymElem(
#          "D3",
#          cen=cell * Vec(0.625, 0.625, 0.625),
#          axis=Vec(1, 1, 1),
#          axis2=Vec(1, -1, 0),
#          col=(0, 1, 0),
#       ),
#       SymElem(
#          "D2",
#          cen=cell * Vec(0.625, 0.500, 0.750),
#          axis=Vec(1, 0, 0),
#          axis2=Vec(0, -1, 1),
#          col=(0, 1, 1),
#       ),
#    ]
#    test_xtal(G, cell, tag="I4132", origin=cell * Vec(0.5, 0.5, 0.5), showshape=0, **kw)
#    # cube( cell*Vec(-0.5,-0.5,-0.5), cell*Vec(0.5,0.5,0.5) )
#    cube(cell * Vec(-0, -0, -0), cell * Vec(1, 1, 1))

# def test_icos(depth=3, cell=30, **kw):
#    # dodec
#    # (+-1, +-1, +-1)
#    # (0, +-1/p, +-p)
#    # (+-1/p, +-p, 0)
#    # (+-p, 0, +-1/p)
#    # icos
#    # (0, +-1, +-p)
#    # (+-1, +-p, 0)
#    # (+-p, 0, +-1)
#    p = (1.0 + sqrt(5.0)) / 2.0
#    q = 1.0 / p
#    v33a = (Vec(1, 1, 1).normalized() - Vec(q, p, 0).normalized()).normalized()
#    v33b = (Vec(1, 1, 1).normalized() - Vec(p, 0, q).normalized()).normalized()
#    v33c = (Vec(1, 1, 1).normalized() - Vec(0, q, p).normalized()).normalized()
#    icsang = angle_degrees(Vec(1, 1, 1), V0, v33a)
#    tetang = 180.0 - math.acos(-1.0 / 3.0) * 180.0 / math.pi
#    print(icsang, tetang)
#    delta_deg = icsang - tetang
#    v33a = RAD(v33a.cross(Vec(1, 1, 1)), delta_deg) * v33a
#    v33b = RAD(v33b.cross(Vec(1, 1, 1)), delta_deg) * v33b
#    v33c = RAD(v33c.cross(Vec(1, 1, 1)), delta_deg) * v33c
#    print(angle_degrees(Vec(1, 1, 1), V0, v33a))
#    print(angle_degrees(Vec(1, 1, 1), V0, v33b))
#    print(angle_degrees(Vec(1, 1, 1), V0, v33c))
#    G = [
#       # icos 3folds
#       SymElem("C3", Vec(+1, +1, +1)),
#       SymElem("C3", Vec(+1, +1, -1)),
#       SymElem("C3", Vec(+1, -1, +1)),
#       SymElem("C3", Vec(-1, +1, +1)),
#       SymElem("C3", Vec(0, +q, +p)),
#       SymElem("C3", Vec(0, +q, -p)),
#       SymElem("C3", Vec(+q, +p, 0)),
#       SymElem("C3", Vec(+q, -p, 0)),
#       SymElem("C3", Vec(+p, 0, +q)),
#       SymElem("C3", Vec(-p, 0, +q)),
#       SymElem("C5", Vec(0, +p, 1)),
#       SymElem("C5", Vec(0, -p, 1)),
#       SymElem("C5", Vec(+p, 1, 0)),
#       SymElem("C5", Vec(-p, 1, 0)),
#       SymElem("C5", Vec(1, 0, +p)),
#       SymElem("C5", Vec(1, 0, -p)),
#       # tet
#       SymElem("C3", Vec(1, 1, 1), cen=cell * Vec(1, 1, 1)),
#       SymElem("C3", v33a, cen=cell * Vec(1, 1, 1)),
#       SymElem("C3", v33b, cen=cell * Vec(1, 1, 1)),
#       SymElem("C3", v33c, cen=cell * Vec(1, 1, 1)),
#       SymElem("C2", v33a + v33b, cen=cell * Vec(1, 1, 1)),
#       SymElem("C2", v33b + v33c, cen=cell * Vec(1, 1, 1)),
#       SymElem("C2", v33c + v33a, cen=cell * Vec(1, 1, 1)),
#    ]
#    symtrie = generate_sym_trie(G, depth=depth)
#    # symtrie.visit(print_node)
#    nodes = []
#    buildcgo = BuildCGO(nodes=nodes, origin=cell * Vec(0.5, 0.5, 0.5), showelems=True, **kw)
#    symtrie.visit(buildcgo)
#    buildcgo.show()

# # def test_F23_C2_T(cell=100,**kw):
# #   #delete all ; run ~/pymolscripts/symgen.py; test_F23_C2_T(cell=100, depth=2)
# #   G = [   SymElem( 'C2', cen=cell* Vec ( 0.0,0.0,0.0 ) , axis=Vec(0,0,1), #axis2=Vec(-,-,-)
# #           ),
# #           SymElem( 'T', cen=cell*Vec(0.25,0.25,0.0), axis=Vec(0.57735,-0.57735,0.57735), #axis2=Vec(,,)
# #           ),
# #
# #          ]
# #   test_xtal(G,cell,tag='test_F23_C2_T',**kw)
# #
# #

# # from una
# def test_F23_C2_T(cell=100, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P4( depth=4, cell=100,
#    # symdef_scale=0.000001, generic_names=1 )
#    G = [
#       SymElem("T", cen=cell * Vec(0.0, 0.0, 0.0), axis=Vec(-1, 1, 1)),
#       SymElem("C2", cen=cell * Vec(0.25, 0.25, 0.0), axis=Vec(0, 0, 1)),
#    ]
#    test_xtal(G, cell, tag="test_F23_C2_T", **kw)

# ## ##
# def test_P321_C3_C2(cell=100, **kw):
#    # delete all ; run ~/pymolscripts/symgen2.py; test_P312_D3_D3(cell=100,
#    # depth=2)
#    G = [
#       SymElem(
#          "C3",
#          cen=cell * (Vec(0, 0, 0)),
#          axis=Vec(0, 0, 1),
#          axis2=Vec(0.866025, 0.5, 0),
#       ),
#       SymElem(
#          "C3",
#          cen=cell * Vec(-0.333333, 0.333333, -0.5),
#          axis=Vec(0, 0, 1),
#          axis2=Vec(0.866025, -0.5, 0),
#       ),
#    ]
#    test_xtal(G, cell, tag="test_P321_C3_C3", **kw)

# ## ##

# def test_P321_C3_C2_will(cell=100, **kw):
#    # delete all ; run ~/pymolscripts/symgen2.py; test_P312_D3_D3(cell=100,
#    # depth=2)
#    G = [
#       SymElem(
#          "C3",
#          cen=cell * (Vec(0, 0, 0)),
#          axis=Vec(0, 0, 1),  # axis2=Vec(0.866025,0.5,0)
#       ),
#       # SymElem('C2', cen=cell * Vec(0, 1, 0), axis=Vec(1, 0, 0),  #
#       # axis2=Vec(0.866025,-0.5,0)
#       SymElem(
#          "C2",
#          cen=cell * Vec(0, 1, 0),
#          axis=Vec(0, 0, 1),  # axis2=Vec(0.866025,-0.5,0)
#       ),
#       # SymElem( 'Z', cen=cell*Vec( 0,0,2 ) ),
#    ]
#    test_xtal(G, cell, tag="test_P321_C3_C2", **kw)

# def test_I432_C2_C4(cell=200, **kw):
#    # delete all; run ~/pymol/symgen.py; test_P4( depth=4, cell=100, symdef_scale=0.000001, generic_names=1 )
#    G = [
#       # SymElem("C4", cen=cell * Vec(0, 0, -0.75), axis=Vec(0, 0, 1)),
#       # SymElem("C2", cen=cell * Vec(0, -0.25, -0.5), axis=Vec(-1, 0.0, 1)),
#       # SymElem("C4", cen=cell * Vec(0, 0, 0), axis=Vec(0, 0, 1)),
#       # SymElem("C2", cen=cell * Vec(0.5, 0.25, 0), axis=Vec(1, 0, 1)),
#       SymElem("O", cen=cell * Vec(0, 0, 0), axis=Vec(0, 0, 1)),
#       # SymElem(
#       # "D2", cen=cell * Vec(0, 0.25, 0.5), axis=Vec(1, 0, 1), axis2=Vec(0, 1, 0)
#       # ),
#       SymElem("C2", cen=cell * Vec(0.5, 0.25, 0), axis=Vec(1, 0, 1)),
#       SymElem("C2", cen=cell * Vec(0.5, 0.25, 0), axis=Vec(-1, 0, 1)),
#    ]
#    test_xtal(G, cell, tag="test_I432_24", **kw)
#    cube(Vec(0, 0, 0), cell * Vec(1, 1, -1))

# def test_P4m_D2_C4(cell=1, **kw):
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C4", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P4m_D2_C4", **kw)

# def test_P4m_D4_C2(cell=1, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C2", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P4m_D4_C2", **kw)

# def test_P4m_D4_C4(cell=1, **kw):
#    G = [
#       SymElem("D4", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C4", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P4m_D4_C4", **kw)

# def test_P3m_D3_C3(cell=1, **kw):
#    G = [
#       SymElem("D3", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C3", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P3m_D3_C3", **kw)

# def test_P6m_D2_C3(cell=1, **kw):
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C3", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P6m_D2_C3", **kw)

# def test_P6m_D2_C6(cell=1, **kw):
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C6", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P6m_D2_C6", **kw)

# def test_P6m_D3_C6(cell=1, **kw):
#    G = [
#       SymElem("D3", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C6", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P6m_D3_C6", **kw)

# def test_P6m_D3_C2(cell=1, **kw):
#    G = [
#       SymElem("D3", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("C2", cen=cell * Vec(1, 0.0, 0.0)),
#    ]
#    test_xtal(G, cell, tag="P6m_D3_C2", **kw)

# def test_F222_D2_D2(cell=100, **kw):
#    G = [
#       SymElem("D2", cen=cell * Vec(0.0, 0.0, 0.0)),
#       SymElem("D2", cen=cell * Vec(1, 1, 1)),
#       SymElem("D2", cen=cell * Vec(1, 1, -1)),
#    ]
#    test_xtal(G, cell, tag="F222_D2_D2", **kw)
