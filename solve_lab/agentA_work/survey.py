import sys,os,json; sys.path.insert(0,'.')
import env, lib as L
S10='/home/user/integer_solver/solve_lab/s10'
rows=[]
for f in sorted(os.listdir(S10)):
    if not f.endswith('.json'): continue
    p=os.path.join(S10,f)
    try: d=json.load(open(p))
    except Exception: continue
    if not isinstance(d,dict) or len(d)<1000: continue
    try: v=L.load(p)
    except Exception as e: continue
    av=L.all_atom_values(v); fe=L.failing_eqs(av)
    nz=[a for a in range(L.NA) if av[a]]
    rows.append((L.NEQ-len(fe), f, len(nz), sorted(fe)[:12], nz[:14]))
rows.sort(reverse=True)
for s,f,n,fe,nz in rows:
    print('%6d  %-28s nzatoms=%-4d fail=%s'%(s,f,n,fe))
    print('          nz: %s'%nz)
