"""GLOBAL mod-p forward evaluation.

Every handle enters as (free var)*(wire), and every wire equals p, so mod p the entire
quotient-witness apparatus disappears and the instance becomes a plain circuit over GF(p):

    free inputs -> gates (each gate determines its output mod p) -> checks must vanish mod p.

If a state has every atom == 0 (mod p), the residue of each equation is p*r and the p-quantised
handles absorb it exactly.  So the mod-p layer is where the real problem lives, and this is the
first time it is evaluated globally rather than in a neighbourhood.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
import fw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

def evalp(Pp, v):
    s=0
    for m,c in Pp.items():
        t=c%P
        for u in m:
            t=t*v[u]%P
            if t==0: break
        s=(s+t)%P
    return s

def coefp(a, t, v):
    """d(atom)/dt at v, mod p (t must occur linearly)"""
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc%P
        for u in m:
            if u!=t: term=term*v[u]%P
        c=(c+term)%P
    return c

def solvep(a, t, v):
    c=coefp(a,t,v)
    if not c: return None
    old=v[t]; v[t]=0
    rest=evalp(L.polys[a], v)
    v[t]=old
    return (-rest)*pow(c,-1,P)%P

def forwardp(v, iters=60):
    """v is a full residue vector; gate outputs get overwritten"""
    undet=[]
    for comp in fw.ORDER:
        if len(comp)==1:
            u=comp[0]
            x=solvep(L.definer[u], u, v)
            if x is None: undet.append(u)
            else: v[u]=x
        else:
            for _ in range(iters):
                ch=False
                for u in comp:
                    x=solvep(L.definer[u], u, v)
                    if x is not None and x!=v[u]: v[u]=x; ch=True
                if not ch: break
    return v, undet

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
    raw=load_raw(src)
    v=[x%P for x in raw]
    t0=time.time()
    v,undet=forwardp(v)
    json.dump([int(x) for x in v], open(os.path.join(HERE,"data","gmp1_state.json"),"w"))
    print(f"mod-p forward eval: {time.time()-t0:.0f}s;  gates that do NOT determine their output "
          f"mod p: {len(undet)} of {len(L.definer)}")
    CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
    bad=[a for a in CHK if evalp(L.polys[a], v)]
    print(f"CHECKS failing mod p: {len(bad)} of {len(CHK)}")
    print("  ", bad[:40])
    badg=[u for u in L.definer if evalp(L.polys[L.definer[u]], v)]
    print(f"gate atoms nonzero mod p after forward eval: {len(badg)} {badg[:10]}")
    # how many equations would fail?
    AVp=[evalp(L.polys[a],v) for a in range(L.NA)]
    F=[e for e in range(L.NEQ)
       if sum(c*AVp[a] for a,c in L.eq_atoms[e][2].items())%P]
    print(f"equations nonzero mod p: {len(F)}  (score ceiling {L.NEQ-len(F)})")
    json.dump({'undet':undet,'bad':bad}, open(os.path.join(HERE,'data','gmp1.json'),'w'))
