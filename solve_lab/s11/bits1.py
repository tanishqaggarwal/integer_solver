"""Everything solid about the 256 real message bits: identity, tree membership, pins."""
import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BOOL=[]
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BOOL.append(ks[0][0])
BOOL=sorted(set(BOOL))
print(f"boolean-checked variables: {len(BOOL)}   free among them: {sum(1 for u in BOOL if u in FREE)}")
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
print(f"REAL (non-inert) bits: {len(real)}   ON: {sorted(u for u in real if base[u]%P==1)}")
# tree membership: which of the four OR-tree leaves does each bit drive?
LEAF=[8599,21839,7304,25956]
grp=collections.defaultdict(list)
for b in real:
    v=list(base); v[b]=(1-base[b])%P; forwardp(v)
    ch=tuple(u for u in LEAF if v[u]!=base[u])
    grp[ch].append(b)
print("\nOR-tree membership (which leaf a bit can drive):")
for k,vv in sorted(grp.items(), key=lambda z:-len(z[1])):
    nm={8599:'A(x8599)',21839:'B(x21839)',7304:'C(x7304)',25956:'D(x25956)'}
    print(f"   {[nm[x] for x in k] if k else 'none (drives no leaf from this state)'}: {len(vv)} bits")
# pins per bit
PIN=collections.defaultdict(list)
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=3: continue
    bb=[m[0] for m,c in Pp.items() if len(m)==1 and m[0] in set(BOOL) and abs(c)>10**60]
    if not bb: continue
    b=bb[0]
    q=[m for m in Pp if len(m)==2 and b in m]
    if not q: continue
    x=[t for t in q[0] if t!=b][0] if q[0][0]!=q[0][1] else q[0][0]
    C=-[c for m,c in Pp.items() if m==(b,)][0]
    PIN[b].append((a,x,C))
npin=collections.Counter(len(PIN[b]) for b in real)
print(f"\nload pins per real bit: {dict(sorted(npin.items()))}")
print(f"total distinct pins on real bits: {sum(len(PIN[b]) for b in real)}")
sz=collections.Counter(len(bin(abs(C))[2:]) for b in real for a,x,C in PIN[b])
print(f"loaded-constant bit-lengths: {dict(sorted(sz.items()))}")
json.dump({'real':sorted(real),'on':sorted(u for u in real if base[u]%P==1),
           'trees':{','.join(map(str,k)):sorted(v) for k,v in grp.items()},
           'pins':{str(b):[[a,x,str(C)] for a,x,C in PIN[b]] for b in real}},
          open(os.path.join(HERE,'data','bits1.json'),'w'))
print("saved data/bits1.json")
