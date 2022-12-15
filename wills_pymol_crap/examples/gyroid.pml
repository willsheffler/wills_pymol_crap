delete all; 
run /home/sheffler/pymol3/symgen_classes.py; run /home/sheffler/pymol3/pymol_util.py; 
run /home/sheffler/pymol3/symgen_classes.py; run /home/sheffler/pymol3/symgen.py;
run /home/sheffler/pymol3/symgen_test.py; 
test_I4132(depth=4,shownodes=0,cell=100,maxrad=80); 
# test_P213(depth=8,shownodes=0,cell=100,maxrad=100); 
run /home/sheffler/pymol3/misc/G222.py; gyroid(100,r=120)
set light, [ -0.3, -0.30, 0.8 ]

zoom

