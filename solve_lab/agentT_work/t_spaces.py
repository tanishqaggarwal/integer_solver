#!/usr/bin/env python3
"""AUDIT T9 -- are agent Q and agent S talking about the SAME 256 booleans?
Q: 256 leaf selectors, fold = group sum, every k in [1,2^256-1] reachable  -> feasible.
S: 256 'cluster booleans', map saturates to one-hot, a20215 != 0 always    -> that target unreachable.
Both say '256'.  Compare the actual variable indices."""
import sys,os,json,collections
Q='/home/user/integer_solver/solve_lab/agentQ_work'
S='/home/user/integer_solver/solve_lab/agentS_work'
lad=json.load(open(os.path.join(Q,'ladder.json')))['ladder']      # exponent -> selector var
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))               # selector var -> [x, y, ...]
Qsel=set(int(v) for v in lad.values())
Qall=set(int(k) for k in qleaf)
os.chdir(S); sys.path.insert(0,S)
import common as C
import engine as E
BOOLS=set(f for f in C.cluster_cone() if C.isbool(f))
cone=set(C.cluster_cone())
print('Q ladder selectors (exponents 0..255, 3 missing) : %d'%len(Qsel))
print('Q qleaf selectors (all 256 decoded leaves)       : %d'%len(Qall))
print('S cluster-cone booleans (bfs.py MOVES minus switch): %d'%len(BOOLS))
print('S full cluster cone (all free vars)              : %d'%len(cone))
print()
print('|Q_ladder  &  S_bools| = %d'%len(Qsel&BOOLS))
print('|Q_qleaf   &  S_bools| = %d'%len(Qall&BOOLS))
print('|Q_qleaf   &  S_cone | = %d'%len(Qall&cone))
print()
print('sample Q ladder selectors :',sorted(Qsel)[:12])
print('sample S cluster booleans :',sorted(BOOLS)[:12])
print()
for a in C.ROWS:
    _,fv=E.cone(a)
    fv=set(fv)
    print('  row a%-6d cone free vars %5d   of which Q leaf selectors: %d'%(a,len(fv),len(fv&Qall)))
