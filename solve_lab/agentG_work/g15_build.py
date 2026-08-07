"""Construct a state from the exact mod-p reduced solve, then lift to Z via handles."""
import os, sys, json, pickle, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gred
from gsym import *
import tools as T

SRC = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
FLIP = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 and sys.argv[2]!='-' else []
OUT  = sys.argv[3] if len(sys.argv)>3 else 'G_build.json'
SYMS = json.load(open('closed_nonbool.json'))

v = L.load(SRC); ad.fwd(v, rounds=6)
for b in FLIP: v[b] = 1 - v[b]
ad.fwd(v, rounds=8)
print('after flip %s: score %d' % (FLIP, L.NEQ-len(L.failing_eqs(L.all_atom_values(v)))), flush=True)

r = gred.reduce_state(v, SYMS)
print('reduced: rank=%d nfree=%d nzc=%d ninc=%d nlin=%d nnon=%d' %
      (r['rank'],r['nfree'],len(r['nzc']),r['ninc'],r['nlin'],r['nnon']), flush=True)
print('residual:', [(a,(g%P if isinstance(g,int) else 'POLY%d'%len(g))) for a,g in r['res']], flush=True)
if r['nzc']: print('nonzero-const checks:', r['nzc'][:10])

M,piv,free = r['M'],r['piv'],r['free']
n=len(SYMS)
t = {c: v[SYMS[c]] % P for c in free}
newv = dict(t)
for row,c in enumerate(piv):
    val = M[row][n] % P
    for c2 in free:
        if M[row][c2]%P: val = (val - M[row][c2]*t[c2]) % P
    newv[c] = val
for c,val in newv.items():
    v[SYMS[c]] = val % P
ad.fwd(v, rounds=10)
av = L.all_atom_values(v)
nzmodp = [a for a in gsym.check_atoms() if av[a] % P]
print('checks nonzero mod p after assignment:', len(nzmodp), nzmodp[:10], flush=True)
nzZ = [a for a in range(L.NA) if av[a]]
print('atoms nonzero over Z:', len(nzZ), 'failing eqs:', len(L.failing_eqs(av)), flush=True)

# ---- integer lift: zero every nonzero check over Z through a solo handle ----
solo = collections.defaultdict(list)
for u in range(L.NVARS):
    if u in L.definer: continue
    ats = L.var_atoms[u]
    if len(ats)==1: solo[ats[0]].append(u)
for rnd in range(8):
    av = L.all_atom_values(v)
    bad = [a for a in range(L.NA) if av[a] and a not in L.atom_out]
    fixed=0
    for a in bad:
        done=False
        for h in solo.get(a,()):
            nv = L.solve_for(a,h,v)
            if nv is not None:
                v[h]=nv; done=True; fixed+=1; break
        if not done:
            for u in sorted(L.avars[a]):
                if u in L.definer: continue
                nv=L.solve_for(a,u,v)
                if nv is not None and nv!=v[u]:
                    pass  # only solo handles: other free inputs would perturb other atoms
    ad.fwd(v, rounds=4)
    av=L.all_atom_values(v)
    bad2=[a for a in range(L.NA) if av[a]]
    print('  lift round %d: fixed %d, nonzero atoms %d, failing %d' % (rnd,fixed,len(bad2),len(L.failing_eqs(av))), flush=True)
    if not fixed: break
av=L.all_atom_values(v)
fail=L.failing_eqs(av)
print('FINAL score %d ; nonzero atoms %d' % (L.NEQ-len(fail), len([a for a in range(L.NA) if av[a]])))
print('failing:', fail[:30])
L.save(v, os.path.join('/home/user/integer_solver/solve_lab/agentG_work',OUT))
print('saved', OUT)
