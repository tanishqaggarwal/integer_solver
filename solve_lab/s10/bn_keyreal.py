"""bn_keyreal: when a key bit moves off {0,1}, do NON-boolean atoms become nonzero?

If yes, that bit can never take part in a pure boolean-carrier configuration
(the cone model assumes every non-boolean atom vanishes), so the mixed-regime
cone is unrealisable through it regardless of the LP verdict.
"""
import os, sys, json
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map(); BA=set(bools)
BASE={22229,22230,35758,35759,35760,35761,35762}
v0=L.load(B.BEST); BG=B.broken_gates(v0)
KEY=[2081,4287,11368,13195]
out={}
for u in KEY:
    for x in (2,-1,3):
        w=list(v0); w[u]=x
        B.fwdb(w,BG,1)
        s,fl,av=B.score(w)
        nz=[a for a in range(L.NA) if av[a]]
        newnz=[a for a in nz if a not in BASE]
        nb=[a for a in newnz if a not in BA]
        bb=[a for a in newnz if a in BA]
        print(f'x_{u}={x}: score {s}  new nonzero atoms {len(newnz)} '
              f'(boolean {len(bb)}, NON-boolean {len(nb)})',flush=True)
        print(f'    boolean: {bb[:12]}')
        print(f'    NON-boolean: {nb[:12]}')
        out[f'{u}={x}']={'score':s,'bool':bb,'nonbool':nb}
json.dump(out, open(os.path.join(HERE,'bn_keyreal.json'),'w'))
print('saved bn_keyreal.json')
