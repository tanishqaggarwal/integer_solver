"""For every free leaf bit: its pin atoms, their handle-definition atoms, and the union of
equations that breaking those handles would touch."""
import sys, json, re, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close3 import *
P=2**256-2**32-977
W='/home/user/integer_solver/solve_lab/agentC_work/'
lp={int(k):v for k,v in json.load(open(W+'leafpts2.json')).items()}
v=[0]*L.NVARS; forward(v)
wires={u for u in range(L.NVARS) if v[u]==P}
NUM=re.compile(r'^x_(\d+) \* \(x_(\d+) - (-?\d+)\)(?: - (\d+) \*)? (?:\* )?')
rows=[]
for b in sorted(lp):
    pins=[]
    for a in L.var_atoms[b]:
        if a in L.atom_out: continue
        s=L.atom_src[a]
        m=re.match(r'^x_%d \* \(x_(\d+) - (-?\d+)\) - (?:(\d+) \* )?x_(\d+)$'%b, s)
        if m:
            pins.append(dict(atom=a,X=int(m.group(1)),C=int(m.group(2)),mult=int(m.group(3) or 1),H=int(m.group(4))))
    if len(pins)<2: continue
    eqs=set(); hd=[]
    ok=True
    for pn in pins:
        d=outs.get(pn['H'])
        if d is None: ok=False; break
        hd.append((d,len(L.atom2eq.get(d,{}))))
        eqs|=set(L.atom2eq.get(d,{}))
    if not ok: continue
    rows.append((len(eqs),sum(x[1] for x in hd),b,[x[0] for x in hd],[x[1] for x in hd],len(pins)))
rows.sort()
print('bits with parsed pins:',len(rows))
print('union-of-equations distribution:',collections.Counter(r[0] for r in rows))
for r in rows[:20]:
    print('  union=%-3d sum=%-3d bit=x_%-6d handles=%s eqs=%s npins=%d'%r)
json.dump([[r[0],r[2],r[3],r[4]] for r in rows],open(W+'bitcost.json','w'))
