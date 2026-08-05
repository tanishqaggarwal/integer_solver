#!/usr/bin/env python3
"""Independent exact verifier: re-evaluates every equation in EQUATIONS.txt over Z.
Usage: python3 verify.py FILE.json    (FILE may hold free inputs or a full assignment)"""
import json, re, sys, time
path=sys.argv[1] if len(sys.argv)>1 else 'SOLUTION.json'
d=json.load(open(path))
NV=38748
vals=[0]*NV
for k,v in d.items():
    i=int(k[2:]) if k.startswith('x_') else int(k)
    vals[i]=int(v)
# if it's a free-input file, complete it via the harness forward pass
if len(d) < NV:
    sys.path.insert(0,'/home/user/integer_solver/solve_lab')
    import os; cwd=os.getcwd(); os.chdir('/home/user/integer_solver/solve_lab')
    import heal_harness as H
    for v in H.freeinp: H.val[v]=vals[v]
    H.forward()
    vals=[H.val[i] for i in range(NV)]
    os.chdir(cwd)
    print(f'(completed {len(d)} free inputs -> full assignment via forward pass)')
pat=re.compile(r'x_(\d+)')
ok=0; bad=[]; t0=time.time(); n=0
with open('/home/user/integer_solver/EQUATIONS.txt') as f:
    for i,line in enumerate(f):
        line=line.strip()
        if not line: continue
        n=i+1
        if line.endswith('= 0'): line=line[:-3].rstrip()
        r=eval(pat.sub(lambda m:'v['+m.group(1)+']', line),{'v':vals,'__builtins__':{}})
        if r==0: ok+=1
        elif len(bad)<15: bad.append(i)
print(f'satisfied {ok}/{n}   ({time.time()-t0:.0f}s)')
if bad: print('failing equations:',bad)
