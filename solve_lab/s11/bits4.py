import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
P=L.P; sys.set_int_max_str_digits(400000)
D=json.load(open(os.path.join(HERE,'data','bits1.json')))
PIN={int(k):v for k,v in D['pins'].items()}
real=sorted(PIN)
pat=collections.Counter()
for b in real:
    ms=[]
    for a,x,C in PIN[b]:
        Pp=L.polys[a]
        z=[abs(c) for m,c in Pp.items() if len(m)==1 and m[0]!=b]
        ms.append(z[0] if z else None)
    pat[tuple(sorted(1 if m==1 else 0 for m in ms))]+=1
print("per-bit pin multiplier pattern (1 = m==1):", dict(pat))
# handle behind each pin
hk=collections.Counter()
for b in real[:60]:
    for a,x,C in PIN[b]:
        Pp=L.polys[a]
        hs=[m[0] for m,c in Pp.items() if len(m)==1 and m[0]!=b]
        if not hs: continue
        h=hs[0]; d=L.definer.get(h)
        if d is None: hk['free']+=1; continue
        mono=[m for m in L.polys[d] if len(m)==2]
        hk['wire*free' if mono else 'other']+=1
print("what sits behind the pin handle (sample of 120 pins):", dict(hk))
flips=json.load(open(os.path.join(HERE,'data','gmp16.json')))
print("\nsingle-flip landscape from the 4-failing base (failing checks mod p):")
print("  ", dict(sorted(collections.Counter(f for f,b,o in flips).items())))
print("  best flips:", [(b,f) for f,b,o in sorted(flips)[:6]])
tr=json.load(open(os.path.join(HERE,'data','bits_trees.json')))
for k,v in tr.items(): print(f"  tree {k}: {len(v)} bits")
