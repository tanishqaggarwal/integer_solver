"""Rigorous: with S = FAILS + {9123,18673} and EVERY variable whose footprint lies inside S
as a free knob, exhaustively test all 5-subsets of S for integer solvability."""
import sys, itertools, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from search2 import knobs_of, build_forms, FAILS
from dioph2 import solve_int
S=sorted(set(FAILS)|{9123,18673})
V=knobs_of(S)
print('S =',S)
print('free knobs (footprint fully inside S):',V)
forms=build_forms(S,V)
Vl=sorted(V)
for k in (6,5):
    ok=0; tot=0
    for sub in itertools.combinations(S,k):
        tot+=1
        M=[[forms[e][1].get(v,0) for v in Vl] for e in sub]
        r=[-forms[e][0] for e in sub]
        if solve_int(M,r) is not None:
            ok+=1; print('SOLVABLE',k,sub)
    print(f'size {k}: {ok}/{tot} integer-solvable')
