import sys, os, json, collections, pickle
F='/home/user/integer_solver/solve_lab/agentV_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
H=pickle.load(open('handles.pkl','rb')); appearP=H['appearP']
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
for w in (36193,35256,10424,27436,20930,30632,19326,28825):
    print('x%d in %d atoms:'%(w,len(appearP.get(w,[]))))
    for a in appearP.get(w,[]): print('     ',a[:240])
