#!/usr/bin/env python3
import json,pickle,re
from propagate import NVARS
p=2**256-2**32-977
D=pickle.load(open('wire_data.pkl','rb')); wire=D['wire']
BP=pickle.load(open('bitpolys.pkl','rb')); polys=BP['polys']; BE=BP['BE']; BITS=BP['BITS']
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v=vA[:]
for y,s in wire.items(): v[y]=s*1
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'v':v,'__builtins__':{}}
bitvals={b:vA[b] for b in BITS}
def poly_eval(P):
    val=P['c0']
    for b,l in P['lin'].items(): val+=l*bitvals[b]
    for b,q in P['quad'].items(): val+=q*bitvals[b]*bitvals[b]
    for (a,bb),cr in P['cross'].items(): val+=cr*bitvals[a]*bitvals[bb]
    return val
bad=0
for i in BE:
    actual=eval(eqcode[i],ns)
    recon=poly_eval(polys[i])
    if actual!=recon:
        bad+=1
        if bad<=6: print(f"MISMATCH eq {i}: actual%p={actual%p}, recon%p={recon%p}")
print(f"extraction check: {len(BE)-bad}/{len(BE)} match exactly at bits=agentA")
# also check the inconsistent eqs contain what vars beyond bits+wire
INC=[32895,33012,33633,26785,32300,36767,37257,37666]
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
for i in INC:
    vs=set(int(m) for m in VAR.findall(lines[i]))
    nonbit_free=[x for x in vs if x in freeinp and x not in BITS]
    inwire=[x for x in vs if x in wire]
    print(f"  incons eq {i}: #vars={len(vs)}, nonbit-free={len(nonbit_free)}, wire-members={len(inwire)}, bits={sorted(vs&set(BITS))}")
