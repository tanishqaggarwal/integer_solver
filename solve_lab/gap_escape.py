import heal_harness as H
import json, glob, pickle
p=H.p
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
def av(i):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
gapatoms=[20862,20864,42669]
print("=== gap atoms 20862(G1),20864(G2),42669 across saved solutions (mod p) ===")
for fn in sorted(glob.glob('best/*.json')):
    try: d=H.loadd(fn)
    except: continue
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    vals=[av(i) for i in gapatoms]
    print(f"{fn.split('/')[-1]:40s} nfail={len(H.fails()):4d}  G1={'0' if vals[0]==0 else 'NZ'}  G2={'0' if vals[1]==0 else 'NZ'}  42669={'0' if vals[2]==0 else 'NZ'}")

# Now cert linearity at 39022
d=H.loadd('best/new_instance_partial_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
from collections import defaultdict
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def incr(w,nv):
    H.val[w]=nv
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],{'v':H.val,'__builtins__':{}})
cd=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/cert22.pkl','rb'))
cert=cd['cert']
def free_anc_atom(i):
    s=set()
    for m,c in ATOMS[i]:
        for v in m:
            if v in H.freeinp: s.add(v)
            else: s|=H.anc.get(v,set())
    return s
nlin=0; lin=0; nonlin=[]
for a,mv in cert:
    fa=sorted(free_anc_atom(a))
    isl=True
    for f in fa[:60]:
        b0=av(a); incr(f,base[f]+1); b1=av(a); incr(f,base[f]+2); b2=av(a); incr(f,base[f])
        if (b1-b0)%p!=(b2-b1)%p: isl=False; break
    if isl: lin+=1
    else: nonlin.append(a); nlin+=1
print(f"\n39022 cert ({len(cert)} atoms): affine={lin}, nonlinear={nlin}")
print("nonlinear cert atoms:",nonlin)
for a in nonlin: print(f"   {a}: {reprs[a][:90]}")
