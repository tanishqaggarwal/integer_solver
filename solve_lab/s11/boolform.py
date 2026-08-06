import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
HERE=os.path.dirname(os.path.abspath(__file__))
LD=json.load(open(os.path.join(HERE,'data','loads.json')))['loads']
BITS=set(int(b) for b in LD)
def fmt(a, lim=170):
    parts=[]
    for mm,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not mm) else ('+' if c==1 else '-')
        if mm: s+='*'.join('x%d'%u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]
# find the SMALLEST atoms containing bit^2
cand=[]
for a in range(L.NA):
    Pp=L.polys[a]
    sq=[mm[0] for mm,c in Pp.items() if len(mm)==2 and mm[0]==mm[1] and mm[0] in BITS]
    if sq: cand.append((len(Pp), a, sq))
cand.sort()
print("smallest atoms with bit^2:")
for n,a,sq in cand[:6]:
    print(f"  a{a} ({n} monomials, {len(L.atom2eq.get(a,{}))} eqs, {'CHECK' if L.atom_out.get(a) is None else 'GATE'}): {fmt(a)}")
